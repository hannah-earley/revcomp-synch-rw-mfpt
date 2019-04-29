#!/usr/bin/env python3
import os
import os.path
import signal
import sys
import json
import time
import datetime
import fcntl
import socket
import subprocess
import contextlib
import argparse
import traceback
import multiprocessing
import pathlib
import email.message
from functools import wraps
import importlib.util

try:
    import psutil
    psutil_ = True
except ImportError:
    psutil_ = False

EXT_JOB = '.job'
EXT_STAT = '.status'

def no_zombies(f):
    @wraps(f)
    def wrapped(*args, **kwargs):
        os.setpgrp()
        try:
            return f(*args, **kwargs)
        except:
            traceback.print_exc()
            os.killpg(0, signal.SIGKILL)
    return wrapped

def handler_help(parser):
    def handler(args):
        parser.print_help()
        sys.exit(1)
    return handler

def handler_none(args):
    raise NotImplementedError

def handler_enq(args):
    qdir = args.q
    jset = args.jobset
    tmpl = args.template
    jdir = os.path.join(qdir, jset)
    jtpl = os.path.join(jdir, tmpl)

    spec = importlib.util.spec_from_file_location("template", jtpl)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)

    jobbs = set()

    for jobn,jobd in mod.generate():
        print(jobn)
        jobf = jobn + EXT_JOB
        jobp = os.path.join(jdir, jobf)
        jobj = Job(jobn, jobp)
        save = lambda f: json.dump(jobd, f, sort_keys=True, indent=4)

        jobb = jobj.filebase(jobd['options'])
        if jobb in jobbs:
            print("  Job would refer to duplicate data file!: %s" % jobb)
            continue
        jobbs.add(jobb)

        try:
            with open(jobp, 'x') as f:
                save(f)
        except FileExistsError:
            try:
                with jobj.lock() as f:
                    f.truncate()
                    save(f)
            except Job.LockedOut:
                with open(jobp, 'r') as f:
                    if json.load(f) != jobj:
                        print("  Can't update -- job being executed!")

def get_stats_(jsdir, base):
    outp = ""
    sep = '  '

    for root, _, jobs in os.walk(jsdir):
        path = pathlib.PurePath(root).relative_to(jsdir)
        parts = path.parts
        lvl = len(parts)
        jobs.sort()

        if lvl:
            if parts[-1] == '__pycache__':
                continue
            outp += sep*lvl + parts[-1] + ':' + "\n"
        else:
            outp += base + ':' + "\n"
        for job in jobs:
            jobn, ext = os.path.splitext(job)
            jobp = os.path.join(root, job)
            if ext != EXT_JOB:
                continue

            jobj = Job(job, jobp)
            outp += sep*(lvl+1) + jobn + ' : ' + jobj.get_status() + "\n"

    return outp

def get_stats(args):
    qdir = args.q
    sets = args.jobset
    if not sets:
        return get_stats_(qdir, 'queue')
    else:
        outp = ""
        for set_ in sets:
            outp += get_stats_(os.path.join(qdir, set_), set_)
        return outp


def handler_stat(args):
    eml = args.email
    int_ = args.interval
    subj = 'Re: RWX Batch Update'
    while True:
        stats = get_stats(args)

        print('[%s] ' % datetime.datetime.now(), end='')
        if eml is None:
            print()
            print(stats)
        else:
            print('emailing to', eml, '...')
            msg = email.message.EmailMessage()
            msg['Subject'] = subj
            msg['From'] = 'RWX Batch <wje25@maths.cam.ac.uk>'
            msg['To'] = eml
            msg.set_content(stats)
            p = subprocess.Popen(["/usr/sbin/sendmail", "-t", "-oi"], stdin=subprocess.PIPE)
            p.communicate(msg.as_string().encode())

        if not int_:
            break
        time.sleep(int_)

class ProcInfo:
    def __init__(self, pid, perc=True):
        self.pid = pid
        self._cpus = multiprocessing.cpu_count()
        self.working = False
        self.perc = perc

        if psutil_:
            self.proc = psutil.Process(pid)
            self.procs = {}
            self.update_procs()
            if self.perc:
                for p in self.procs.values():
                    p.cpu_percent()
                time.sleep(1)
            self.working = True

    def update_procs(self):
        if not psutil_:
            return
        procsn = {}
        try:
            for p in [self.proc] + self.proc.children(recursive=True):
                if p.status().lower() != 'running':
                    continue
                pid = p.pid
                if pid in self.procs:
                    procsn[pid] = self.procs[pid]
                else:
                    procsn[pid] = p
        except psutil.NoSuchProcess:
            pass
        self.procs = procsn


    def cpup(self, update=True):
        if not self.working or not self.perc:
            return
        if update:
            self.update_procs()
        return sum(p.cpu_percent() for p in self.procs.values())

    def cpus(self, update=True):
        if not self.working:
            return
        if update:
            self.update_procs()
        try:
            aff = set()
            for p in self.procs.values():
                aff.update(p.cpu_affinity())
            return len(aff)
        except AttributeError:
            return self._cpus

    def cpusp(self, update=True):
        if self.working and update:
            self.update_procs()
        return self.cpus(False), self.cpup(False)

def cpus():
    return ProcInfo(None,False).cpus()

class Job:
    class TryAgain(Exception): pass
    class LockedOut(Exception): pass

    def __init__(self, job, path):
        self.loaded = False
        self.name = job
        self.path = path
        self.path_stat = os.path.splitext(path)[0] + EXT_STAT

    def load(self, force=False):
        if self.loaded and not force:
            return

        with open(self.path, 'r') as f:
            descriptor = json.load(f)

        self.opts = descriptor['options']
        self.target = descriptor['target']
        self.reqs = descriptor['requirements']
        self.skip = 0

        if isinstance(self.opts, str):
            self.opts = self.opts.split()
        if 'skip' in self.target:
            self.skip = self.target['skip']

        self.filebase = self.filebase(self.opts)
        self.loaded = True

    def needload(f):
        @wraps(f)
        def wrapped(self, *args, **kwargs):
            self.load()
            return f(self, *args, **kwargs)
        return wrapped

    @classmethod
    def filebase(cls, opts):
        opts = ['./job.sh'] + opts + ['-P', '-i', '0']
        outp = subprocess.Popen(opts, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        sout, serr = outp.communicate()
        lines = sout.decode().strip().split('\n')
        assert len(lines) == 1, "Job.filebase: Invalid output from ./job.sh"
        base = lines[0].strip()
        assert base, "Job.filebase: No file base found"
        return base

    @needload
    def output_count(self):
        path = self.filebase + ".csv"
        n = -self.skip
        try:
            with open(path, 'r') as f:
                for line in f:
                    if line.strip()[0] == '#':
                        continue
                    n += 1
        except FileNotFoundError:
            pass
        return n

    @needload
    def get_error(self):
        import refine
        path = self.filebase + ".csv"
        try:
            with open(path, 'r') as f:
                _, err = refine.refine(f, self.skip)
                return err
        except (FileNotFoundError, ZeroDivisionError):
            return float('inf')

    @needload
    def get_times(self):
        import rawfine
        pin = self.filebase + ".log"
        pdat = self.filebase + ".dat"
        try:
            for rec in rawfine.read_log(pin, pdat): pass
            _, est_time = rec['wall_']
            prog = rec['prog']
            return prog, est_time
        except:
            return '??%', '??'

    @needload
    def task_size(self):
        curr = self.output_count()
        if 'count' in self.target:
            count = self.target['count']
            return max(count - curr, 0)

        elif 'chunk' in self.target:
            chunk = self.target['chunk']
            min_ = self.target.get('min', 0)
            max_ = self.target.get('max', float('inf'))
            prec = self.target['prec']
            err = self.get_error()

            if (curr >= min_ and err < prec) or curr >= max_:
                return 0
            return max(min(chunk, max_-curr), 0)

    @needload
    def can_run(self):
        if 'cpu' in self.reqs:
            min_cpu = self.reqs['cpu']
            cur_cpu = cpus()
            if cur_cpu < min_cpu:
                print(self.path, ": Too few cpus!")
                return False
        return True

    @contextlib.contextmanager
    def lock(self):
        with open(self.path, 'r+') as fh:
            fd = fh.fileno()
            try:
                fcntl.lockf(fd, fcntl.LOCK_EX | fcntl.LOCK_NB)
            except OSError:
                raise Job.LockedOut()

            try:
                yield fh
            finally:
                fcntl.lockf(fd, fcntl.LOCK_UN)

    @needload
    def status(self):
        host = socket.gethostname()
        if len(host) > 25:
            host = host[:22] + '...'

        cpuu = ''
        n = self.output_count()

        try:
            pi = self.proc_info
        except AttributeError:
            pass
        else:
            cpus, cpup = pi.cpusp()
            if cpus is not None and cpup is not None:
                cpuu = ' cpu:%.2f/%d' % (cpup/100.0,cpus)
                try:
                    cpuu += ' (%.1f%%)' % (cpup/cpus)
                except ZeroDivisionError:
                    pass

        stat = 'outs:%d' % n
        if 'count' in self.target:
            stat += '/%d' % self.target['count']
        elif 'chunk' in self.target:
            prec = self.target['prec']
            if 'max' in self.target:
                stat += '/%d' % self.target['max']

        err = self.get_error()
        stat += ' err:%g' % err
        if 'prec' in self.target:
            stat += '/%g' % self.target['prec']

        prog, wall = self.get_times()
        stat += ' prog:%s ~time:%s' % (prog,wall)

        now = datetime.datetime.now()
        ts = now.strftime('%y-%m-%d %H:%M')
        return '%s [%s]%s %s' % (ts,host,cpuu,stat)

    def record_status(self, msg='RUNNING', disp=False):
        s = self.status()
        _, s_ = s.split('] ')
        with open(self.path_stat, 'w') as f:
            f.write('%s %s' % (msg, s))
        if disp:
            print('%s: %s' % (self.path, s_))

    def get_status(self):
        status = "QUEUED"
        cached = ""
        try:
            with open(self.path_stat) as f:
                cached = f.read()
        except FileNotFoundError:
            pass

        try:
            with self.lock():
                if cached:
                    status = "HALTED"
        except Job.LockedOut:
            status = "ACTIVE"

        return '[%s] %s' % (status, cached)


    @needload
    def execute(self):
        if not self.can_run():
            return
        with self.lock():
            started = False
            while True:
                todo = self.task_size()
                if todo <= 0:
                    if started:
                        self.record_status(msg='DONE')
                    return True
                started = True
                self.record_status(msg='INITIALISING')

                args = ["./job.sh"] + self.opts + ["-i", str(todo)]
                dn = subprocess.DEVNULL
                self.proc = subprocess.Popen(args, stdout=dn, stderr=dn)
                self.proc_info = ProcInfo(self.proc.pid)
                i = 0
                while True:
                    ret = self.proc.poll()
                    if ret is not None:
                        if ret > 0:
                            return
                        break
                    self.record_status(disp=(i%6==0))
                    time.sleep(10)
                    i += 1
                del self.proc_info

def run_job(job, path):
    jobj = Job(job, path)
    try:
        return jobj.execute()
    except Job.LockedOut:
        return # don't try again...
        # return False # try again...
    except Job.TryAgain:
        return False

@no_zombies
def handler_run(args):
    # first, set highest nice level to be cooperative
    # on a shared server, unless explicitly told it's ok
    if not args.mean:
        os.nice(40)

    qdir = args.q

    sets = args.jobset
    if not sets:
        sets = ['']
    setdirs = [os.path.join(qdir, set_) for set_ in sets]

    again = True
    visited = set()
    while again or args.inf:
        again = False
        for dir_ in setdirs:
            for root, _, jobs in os.walk(dir_):
                jobs.sort()
                for job in jobs:
                    _, ext = os.path.splitext(job)
                    jobpath = os.path.join(root, job)
                    if jobpath in visited:
                        continue
                    if ext == EXT_JOB and run_job(job, jobpath) == False:
                        again = True
                    else:
                        visited.add(jobpath)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="")
    parser.add_argument('-q', default='./queue',
                        help="queue directory")
    subparsers = parser.add_subparsers(title='commands', help='run with -h for more info')
    parser.set_defaults(handler=handler_help(parser))

    parser_enq = subparsers.add_parser('enqeue', aliases=['enq'])
    parser_enq.add_argument('jobset')
    parser_enq.add_argument('--template', default='template.py')
    parser_enq.set_defaults(handler=handler_enq)
    
    parser_stat = subparsers.add_parser('status', aliases=['stat'])
    parser_stat.add_argument('-n', '--interval', default=None, type=float,
                             help="repeat at interval")
    parser_stat.add_argument('-e', '--email', default=None, type=str,
                             nargs='?', const='wje25@cam.ac.uk',
                             help="send stats by email")
    parser_stat.add_argument('jobset', nargs='*')
    parser_stat.set_defaults(handler=handler_stat)
    
    parser_run = subparsers.add_parser('run')
    parser_run.add_argument('--mean', default=False, action='store_true',
                            help="don't run nicely (e.g. hpc allocation)")
    parser_run.add_argument('--inf', default=False, action='store_true',
                            help="run indefinitely")
    parser_run.add_argument('jobset', nargs='*')
    parser_run.set_defaults(handler=handler_run)

    args = parser.parse_args()
    args.handler(args)

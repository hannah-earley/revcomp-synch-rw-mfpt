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

EXT_JOB = '.job'
EXT_STAT = '.status'

def handler_help(parser):
    def handler(args):
        parser.print_help()
        sys.exit(1)
    return handler

def handler_none(args):
    raise NotImplementedError

handler_enq = handler_none

handler_stat = handler_none

def proc_info(pid, sample_time=1.0):
    try:
        import psutil
    except ImportError:
        return None
    else:
        p = psutil.Process(pid)
        cs = [p] + p.children(recursive=True)
        for c in cs:
            c.cpu_percent()
        time.sleep(sample_time)
        cpup = sum(c.cpu_percent() for c in cs)

        try:
            aff = set()
            for c in cs:
                aff += set(c.cpu_affinity())
            cpun = len(aff)
        except AttributeError:
            import multiprocessing
            cpun = multiprocessing.cpu_count()

    return cpup, cpun



class Job:
    class TryAgain(Exception): pass
    class LockedOut(Exception): pass

    def __init__(self, job, path):
        self.name = job
        self.path = path
        self.path_stat = os.path.splitext(path)[0] + EXT_STAT
        self.load()

    def load(self):
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

        self.filebase = self.get_filebase()

    def get_filebase(self):
        opts = ['./job.sh'] + self.opts + ['-P', '-i', '0']
        outp = subprocess.run(opts, capture_output=True, encoding='utf-8')
        lines = outp.stdout.strip().split('\n')
        assert len(lines) == 1, "Job.get_filebase: Invalid output from ./job.sh"
        base = lines[0].strip()
        assert base, "Job.get_filebase: No file base found"
        return base

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

    def get_error(self):
        import refine
        path = self.filebase + ".csv"
        try:
            with open(path, 'r') as f:
                _, err = refine.refine(f, self.skip)
                return err
        except (FileNotFoundError, ZeroDivisionError):
            return float('inf')

    def task_size(self):
        curr = self.output_count()
        if 'count' in self.target:
            count = self.target['count']
            return max(count - curr, 0)

        elif 'chunk' in self.target:
            chunk = self.target['chunk']
            limit = self.target['limit']
            prec = self.target['prec']
            err = self.get_error()

            if err < prec or curr >= limit:
                return 0
            return max(min(chunk, limit-curr), 0)

    def can_run(self):
        # could check job requirements here...
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
                yield
            finally:
                fcntl.lockf(fd, fcntl.LOCK_UN)

    def status(self, pid=None):
        host = socket.gethostname()
        cpuu = ''
        n = self.output_count()

        if pid:
            cpui = proc_info(pid)
            if cpui is not None:
                cpup, cpun = cpui
                cpuu = ' cpu:%.2f/%d (%.1f%%)' % (cpup/100.0,cpun,cpup/cpun)

        stat = 'outs:%d' % n
        if 'count' in self.target:
            stat += '/%d' % self.target['count']
        elif 'chunk' in self.target:
            prec = self.target['prec']
            err = self.get_error()
            stat += '/%d' % self.target['limit']
            stat += ' err:%g/%g' % (err,prec)

        return '%s [%s]%s %s' % (datetime.datetime.now(),host,cpuu,stat)

    def record_status(self, pid=None, msg='RUNNING', disp=False):
        s = self.status(pid)
        _, s_ = s.split('] ')
        with open(self.path_stat, 'w') as f:
            f.write('%s %s' % (msg, s))
        if disp:
            print('%s: %s' % (self.path, s_))

    def execute(self):
        if not self.can_run():
            return
        with self.lock():
            self.record_status(msg='INITIALISING')
            while True:
                todo = self.task_size()
                if todo <= 0:
                    self.record_status(msg='DONE')
                    return True

                args = ["./job.sh"] + self.opts + ["-i", str(todo)]
                dn = subprocess.DEVNULL
                proc = subprocess.Popen(args, stdout=dn, stderr=dn)
                i = 0
                while proc.poll() is None:
                    self.record_status(proc.pid) #, disp=(i%6==0))
                    time.sleep(10)
                    i += 1

def run_job(job, path):
    jobj = Job(job, path)
    try:
        return jobj.execute()
    except Job.LockedOut:
        return # don't try again...
        # return False # try again...
    except Job.TryAgain:
        return False


def handler_run(args):
    # first, set highest nice level to be cooperative
    # on a shared server, unless explicitly told it's ok
    if not args.mean:
        os.nice(40)

    qdir = args.q

    sets = args.jobset
    if not sets:
        sets = ['']
    setdirs = [qdir + '/' + set_ for set_ in sets]

    again = True
    visited = set()
    while again or args.inf:
        again = False
        for dir_ in setdirs:
            for root, _, jobs in os.walk(dir_):
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
    os.setpgrp()
    try:

        parser = argparse.ArgumentParser(description="")
        parser.add_argument('-q', default='./queue',
                            help="queue directory")
        subparsers = parser.add_subparsers(title='commands', help='run with -h for more info')
        parser.set_defaults(handler=handler_help(parser))

        parser_enq = subparsers.add_parser('enqeue')
        parser_enq.set_defaults(handler=handler_enq)
        
        parser_stat = subparsers.add_parser('status')
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

    except:
        traceback.print_exc()
        os.killpg(0, signal.SIGKILL)
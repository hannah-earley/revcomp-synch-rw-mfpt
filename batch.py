#!/usr/bin/env python3
import os
import sys
import json
import time
import fcntl
import subprocess
import contextlib
import argparse

def handler_help(parser):
    def handler(args):
        parser.print_help()
        sys.exit(1)
    return handler

def handler_none(args):
    raise NotImplementedError

handler_enq = handler_none

handler_stat = handler_none


class Job:
    class TryAgain(Exception): pass
    class LockedOut(Exception): pass

    def __init__(self, job, path):
        self.name = job
        self.path = path
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

    def execute(self):
        with self.lock():
            print(self.opts)
            print(self.target)
            print(self.reqs)
            print(self.skip)
            while True:
                todo = self.task_size()
                if todo <= 0:
                    return True
                time.sleep(1)
                return True






# def run_job(job, path):
#     def err(s):
#         print("%s: %s (%s)" % (job,s,path))
#     with open(path, 'r+') as fh:
#         fd = fh.fileno()
#         try: # acquire exclusive lock
#             fcntl.lockf(fd, fcntl.LOCK_EX | fcntl.LOCK_NB)
#         except OSError:
#             return False
#         try:
#             try: # read job descriptor
#                 deets = json.load(fh)
#             except json.JSONDecodeError:
#                 return err("Invalid job descriptor")

#             print(job, path, deets)
#             try:
#                 jobj = Job(job, path, deets)
#             except Exception as e:
#                 return err("Couldn't initialise job [%r]" % e)
#             time.sleep(1)
#             return True


#         finally: # release lock
#             fcntl.lockf(fd, fcntl.LOCK_UN)
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
                    jobpath = os.path.join(root, job)
                    if jobpath in visited:
                        continue
                    if run_job(job, jobpath) == False:
                        again = True
                    else:
                        visited.add(jobpath)


if __name__ == '__main__':
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
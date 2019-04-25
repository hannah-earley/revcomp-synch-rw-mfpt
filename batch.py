#!/usr/bin/env python3
import os
import sys
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
def handler_run(args):
    # first, set highest nice level to be cooperative
    # on a shared server, unless explicitly told it's ok
    if not args.mean:
        os.nice(40)
    handler_none(args)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="")
    subparsers = parser.add_subparsers(title='commands', help='run with -h for more info')
    parser.set_defaults(handler=handler_help(parser))

    parser_enq = subparsers.add_parser('enqeue')
    parser_enq.set_defaults(handler=handler_enq)
    
    parser_stat = subparsers.add_parser('status')
    parser_stat.set_defaults(handler=handler_stat)
    
    parser_run = subparsers.add_parser('run')
    parser_run.add_argument('--mean', default=False, action='store_true',
                        help="don't run nicely (e.g. hpc allocation)")
    parser_run.add_argument('jobset', nargs='*')
    parser_run.set_defaults(handler=handler_run)

    args = parser.parse_args()
    args.handler(args)
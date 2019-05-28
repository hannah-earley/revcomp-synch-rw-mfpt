#!/usr/bin/env python3
import argparse
import shutil
import textwrap
from functools import wraps
import sys
import math

def run_once(message=None, error=True, warn=True):
    def wrapper(f):
        nonlocal message
        if not message:
            message = "Attempt to run %s multiple times!" % f.__name__
        has_run = False
        ret_val = None

        @wraps(f)
        def wrapped(*args, **kwargs):
            nonlocal has_run, ret_val
            if has_run:
                if error:
                    raise RuntimeError(message)
                if warn:
                    print("Warning: " + message)
                return ret_val
            has_run = True
            ret_val = f(*args, **kwargs)
            return ret_val
        return wrapped

    if callable(message):
        f = message
        message = None
        return wrapper(f)
    return wrapper

def handler_help(parser):
    def handler(args):
        parser.print_help()
        parser.exit()
    return handler

def handler_none(args):
    print(args)
    raise NotImplementedError

class _FuncAction(argparse.Action):
    def __init__(self,
                 option_strings,
                 dest=argparse.SUPPRESS,
                 default=argparse.SUPPRESS,
                 help_=None):
        super().__init__(
            option_strings=option_strings,
            dest=dest,
            default=default,
            nargs=0,
            help=help_)

    def __call__(self, parser, namespace, values, option_string=None):
        f = self._fa_func
        args = self._fa_args
        kwargs = self._fa_kwargs

        if self._fa_more:
            kwargs['parser'] = parser
            kwargs['namespace'] = namespace
            kwargs['values'] = values
            kwargs['option_string'] = option_string

        return f(*args, **kwargs)

def FuncAction(f, *args, more=False, **kwargs):
    # more - get action parameters in kwargs
    name = 'FuncAction<%s>' % f.__name__
    dict_ = dict(_fa_func=f,
                 _fa_args=args,
                 _fa_kwargs=kwargs,
                 _fa_more=more)
    return type(name, (_FuncAction,), dict_)

class CommandParser:
    MAX_FIELD_WIDTH = 24
    FIELD_SPACER = 2

    def __init__(self, *args, **kwargs):
        kwargs['add_help'] = False
        parser = argparse.ArgumentParser(*args, **kwargs)

        self.common_args = []
        subclass = type('CommonParser', (CommonParser,),
            dict(common_args=self.common_args))

        self.parser = parser
        self.root = (parser, args, kwargs)
        self.subparsers = parser.add_subparsers(title='commands',
                metavar='command',
                help='Supply with -h for more information',
                parser_class=subclass)
        self.commands = []
        self.handler(handler_help(self))

        self.long_help = False
        parser.format_help = self.format_help
        help_action = FuncAction(self.action_help, more=True)
        help_group = parser.add_mutually_exclusive_group()
        help_group.add_argument('-h', action=help_action)
        help_group.add_argument('--help', action=help_action)

    def handler(self, handler=handler_none):
        self.parser.set_defaults(handler=handler)
        return handler

    def add_command(self, *args, **kwargs):
        parser = self.subparsers.add_parser(*args, **kwargs)
        self.commands.append((parser, args, kwargs))
        return parser

    def add_argument(self, *args, **kwargs):
        self.common_args.append((args, kwargs))

    def action_help(self, **kwargs):
        if kwargs.get('option_string').endswith('help'):
            self.long_help = True

        self.print_help()
        self.exit()

    def format_help(self):
        help_ = ""
        root, rargs, rkwargs = self.root
        cmds = self.commands
        
        help_ += root.format_usage()
        for cmd,_,_ in cmds:
            help_ += cmd.format_usage()

        lines = []
        lines.append('')

        desc = rkwargs.get('description')
        if desc:
            lines.append(desc)
            lines.append('')

        lines.append('optional arguments:')
        lines.append(('  -h', 'show this help message and exit'))
        lines.append(('  --help', 'show this help message and exit (long form)'))
        lines.append('')

        lines.append('commands:')
        for cmd,cargs,ckwargs in cmds:
            name = cargs[0]
            aliases = ckwargs.get('aliases', [])
            names = [name] + aliases
            chelp = ckwargs.get('help', '')
            lines.append(('  ' + ', '.join([name] + aliases), chelp))

        # calculate column sizes...
        cols, rows = shutil.get_terminal_size((80, 25))
        cols = max(cols - 2, 2 * self.MAX_FIELD_WIDTH + self.FIELD_SPACER)
        field1s = [line[0] for line in lines if isinstance(line,tuple)]
        field1_width = max(len(field1) for field1 in field1s)
        field1_width = min(self.MAX_FIELD_WIDTH, field1_width)
        field2_width = cols - field1_width - self.FIELD_SPACER
        indent = ' ' * field1_width
        spacer = ' ' * self.FIELD_SPACER

        for line in lines:
            try:
                field1, field2 = line
                help_ += field1

                if field2:
                    field1n = len(field1)
                    if field1n > field1_width:
                        help_ += "\n" + indent
                    else:
                        help_ += ' ' * (field1_width - field1n)

                    x, *xs = textwrap.wrap(field2, field2_width)
                    help_ += spacer + x + "\n"
                    for x in xs:
                        help_ += indent + spacer + x + "\n"
                else:
                    help_ += "\n"

            except ValueError:
                line = textwrap.fill(line, cols)
                help_ += line + "\n"


        if self.long_help:
            for cmd,cargs,ckwargs in cmds:
                help_ += "\n---\n\n"
                help_ += cmd.format_help()

        return help_

    def parse_args(self, *args, **kwargs):
        return self.parser.parse_args(*args, **kwargs)

    def print_help(self, *args, **kwargs):
        return self.parser.print_help(*args, **kwargs)

    def exit(self, *args, **kwargs):
        return self.parser.exit(*args, **kwargs)

class CommonParser(argparse.ArgumentParser):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.install_common_args()
        self.handler()

    def install_common_args(self):
        for args, kwargs in self.common_args:
            self.add_argument(*args, **kwargs)

    def handler(self, handler=handler_none):
        self.set_defaults(handler=handler)
        return handler

class peekable:
    def __init__(self, xi):
        self._iter = iter(xi)
        self._fin = False
        self._peeked = False
        self._cache = None

    def __iter__(self):
        return self

    def __next__(self):
        if self._fin:
            raise StopIteration()

        if not self._peeked:
            self.peek()

        self._peeked = False
        return self._cache

    def peek(self):
        if self._peeked:
            return self._cache

        try:
            self._cache = next(self._iter)
            self._peeked = True
            return self._cache
        except StopIteration:
            self._fin = True
            raise

    def safe_peek(self):
        try:
            return True, self.peek()
        except StopIteration:
            return False, None

class LCS:
    def __init__(self, xs, ys):
        self.xs = xs
        self.ys = ys
        self.compute()

    @staticmethod
    def mergeSet(xs, ys):
        xi = peekable(xs)
        yi = peekable(ys)

        while True:
            xp, x = xi.safe_peek()
            yp, y = yi.safe_peek()

            if xp:
                if yp:
                    if x < y:
                        yield next(xi)
                    elif x == y:
                        yield next(xi)
                        next(yi)
                    else:
                        yield next(yi)
                else:
                    yield from xi
                    break
            else:
                yield from yi
                break

    @staticmethod
    def unique_diff(zs, z):
        d1,_,_ = z
        if d1 >= 0:
            return zs + [z]
        lz = len(zs)
        for i in range(lz-1,-1,-1):
            d2,_,_ = zs[i]
            if d2 <= 0:
                return zs[:i+1] + [z] + zs[i+1:]
        return [z] + zs

    def compute(self):
        # Hirschberg-esque algorithm
        # uses O(min(m,n)) space rather than O(mn)
        # we also compute the diff simultaneously

        xs = self.xs
        ys = self.ys

        lx = len(xs)
        ly = len(ys)

        if lx < ly:
            lx, ly = ly, lx
            xs, ys = ys, xs

        ud = LCS.unique_diff
        row = lambda i: [[(0,[(-1,k,x) for k,x in enumerate(xs[:i])])]] + [None]*ly
        memo = [[[(0,[(+1,k,y) for k,y in enumerate(ys[:j])])]
                    for j in range(ly+1)],
                row(1)]

        for i,x in enumerate(xs):
            if i > 0:
                memo = [memo[1], row(i+1)]

            for j,y in enumerate(ys):
                if x == y:
                    memo[1][j+1] = [(n+1,zs+[(0,(i,j),x)]) for (n,zs) in memo[0][j]]
                else:
                    zs1 = [(n,ud(zs, (+1,j,y))) for (n,zs) in memo[1][j]]
                    zs2 = [(n,ud(zs, (-1,i,x))) for (n,zs) in memo[0][j+1]]
                    zss = list(LCS.mergeSet(zs1, zs2))
                    lz = max(n for n,_ in zss)
                    memo[1][j+1] = [(n,zs) for (n,zs) in zss if n == lz]

        self.memo = memo[1]

    def common(self):
        seqs = [[z for (d,_,z) in zs if d == 0] for _,zs in self.memo[-1]]
        return seqs

    def common_str(self):
        return [''.join(zs) for _,zs in self.common()]

    def __len__(self):
        return self.memo[-1][-1][0]

    def diff(self):
        return [zs for _,zs in self.memo[-1]]

    def diff_summary(self):
        return [[z for z in zs if z[0] != 0] for zs in self.diff()]


class py3_cmp:
    def __eq__(self, other):
        return self.__cmp__(other) == 0
    def __ne__(self, other):
        return self.__cmp__(other) != 0
    def __gt__(self, other):
        return self.__cmp__(other) > 0
    def __lt__(self, other):
        return self.__cmp__(other) < 0
    def __ge__(self, other):
        return self.__cmp__(other) >= 0
    def __le__(self, other):
        return self.__cmp__(other) <= 0

def cmp(a, b):
    try:
        return a.__cmp__(b)
    except AttributeError:
        try:
            return -b.__cmp__(a)
        except AttributeError:
            return (a>b) - (a<b)

def cmp_float(a, b, tolerance, tol_abs=None):
    if tol_abs is True:
        tol_abs = tolerance

    try:
        if abs(a - b) / a > tolerance:
            return cmp(a,b)
    except ZeroDivisionError:
        if tol_abs:
            if abs(b) > tol_abs:
                return cmp(a,b)
        else:
            return cmp(a, b)
    return 0

def LineIterator(lines, comment_prefix='#', strip=True):
    for line in lines:
        if strip:
            line = line.strip()
        if line.startswith(comment_prefix):
            continue
        yield line

def SkipIterator(it, skip=0):
    for x in it:
        if skip > 0:
            skip -= 1
        else:
            yield x

def ProgressIterator(it, limit=float('inf'), every=None, file=sys.stderr):
    log = lambda *a,**k: print(*a,**k,file=file,end='',flush=True)
    count = 0

    limit = abs(limit)
    if every is None:
        every = limit * 0.05

    if every == float('inf'):
        yield from it
        return

    stat = lambda: str(count)
    if limit < float('inf'):
        digits = math.ceil(math.log10(limit / every))
        places = max(0, digits - 2)
        fmt = '%.' + str(places) + 'f%%'
        stat = lambda: fmt % (100*count/limit)

    log(stat())
    for x in it:
        yield x
        count += 1
        if count % every == 0:
            log('..' + stat())
        if count >= limit:
            break
    print()

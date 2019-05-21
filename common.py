#!/usr/bin/env python3
import argparse
import shutil
import textwrap
from functools import wraps

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

"""
-- longest common subsequence algorithm
-- inefficient - improve by memoising lcs...

lcs :: Ord a => [a] -> [a] -> [[a]]
lcs (x:xs) (y:ys)
  | x == y = map (x:) (lcs xs ys)
  | x /= y = longest $ mergeSet (lcs (x:xs) ys) (lcs xs (y:ys))
lcs _ _ = [[]]

mergeSet :: Ord a => [a] -> [a] -> [a]
mergeSet [] xs = xs
mergeSet xs [] = xs
mergeSet (x:xs) (y:ys) = case compare x y of
    LT -> x : mergeSet xs (y:ys)
    EQ -> x : mergeSet xs ys
    GT -> y : mergeSet (x:xs) ys

longest :: [[a]] -> [[a]]
longest = go 0 []
  where
    go _ yss [] = reverse yss
    go n yss (xs:xss) =
      let l = length xs in
      case compare l n of
        LT -> go n yss xss
        EQ -> go n (xs:yss) xss
        GT -> go l [xs] xss
"""

def longest_subseq_str(xs, ys, cmp=lambda a,b:a==b):
    lx = len(xs)
    ly = len(ys)

    memo = [[None for _ in range(ly+1)] for _ in range(lx+1)]
    for i in range(lx+1):
        memo[i][0] = {""}
    for j in range(ly+1):
        memo[0][j] = {""}

    for i,x in enumerate(xs):
        for j,y in enumerate(ys):
            if cmp(x, y):
                memo[i+1][j+1] = {zs+x for zs in memo[i][j]}
            else:
                zss = memo[i+1][j].union(memo[i][j+1])
                lz = max(map(len, zss))
                memo[i+1][j+1] = {zs for zs in zss if len(zs) == lz}

    return memo[lx][ly]

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

def longest_subseq(xs, ys):
    lx = len(xs)
    ly = len(ys)

    memo = [[None for _ in range(ly+1)] for _ in range(lx+1)]
    for i in range(lx+1):
        memo[i][0] = [[]]
    for j in range(ly+1):
        memo[0][j] = [[]]

    for i,x in enumerate(xs):
        for j,y in enumerate(ys):
            if x == y:
                memo[i+1][j+1] = [zs+[x] for zs in memo[i][j]]
            else:
                zss = list(mergeSet(memo[i+1][j], memo[i][j+1]))
                lz = max(map(len, zss))
                memo[i+1][j+1] = [zs for zs in zss if len(zs) == lz]

    return memo[lx][ly]

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
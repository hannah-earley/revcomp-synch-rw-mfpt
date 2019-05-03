#!/usr/bin/env python3
import argparse
import shutil
import textwrap

def handler_help(parser):
    def handler(args):
        parser.print_help()
        parser.exit()
    return handler

def handler_none(args):
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
#!/usr/bin/env python3
import os
import sys

def nice_int(n):
    parts = []
    while n > 0:
        n, r = divmod(n, 1000)
        if n > 0:
            s = str(1000+r)[-3:]
        else:
            s = str(r)
        parts.append(s)
    return ','.join(reversed(parts))

class Padder:
    def __init__(self, length, s, padstr):
        excess = length - len(s)
        padlen = len(padstr)
        self._str = s
        self._pad = (padstr * ((excess//padlen)+1))[:excess]

class L(Padder):
    def __str__(self):
        return self._str + self._pad
class R(Padder):
    def __str__(self):
        return self._pad + self._str
class C(Padder):
    def __str__(self):
        pad = self._pad
        l = len(pad)//2
        padl, padr = pad[:l], pad[l:]
        return padl + self._str + padr

class Columns:
    class HLine:
        def __init__(self, col, rep):
            self._counts = col._counts
            self._rep = rep
        def __iter__(self):
            for c in self._counts:
                yield str(L(c, '', self._rep))

    def __init__(self, *alignments, separator='  ', padstr=' '):
        self._alignments = alignments
        self._data = []
        self._separator = separator
        self._padstr = padstr
        self._counts = []

    def __call__(self, *fields):
        while len(fields) > len(self._counts):
            self._counts.append(0)
        while len(fields) > len(self._alignments):
            self._alignments.append(L)
        for n,l in enumerate(map(len,fields)):
            self._counts[n] = max(self._counts[n],l)
        self._data.append(fields)

    def hline(self, rep='='):
        self._data.append(Columns.HLine(self, rep))

    def rows(self):
        def row(r):
            for align,count,field in zip(self._alignments, self._counts, r):
                yield align(count, field, self._padstr)
        row_ = lambda r: self._separator.join(map(str, row(r)))
        yield from map(row_, self._data)

    def __str__(self):
        return '\n'.join(map(str, self.rows()))

if __name__ == '__main__':
    try:
        dir_ = sys.argv[1]
        ext = 'csv'
        if len(sys.argv) >= 3:
            ext = sys.argv[2]
        assert(ext in ['csv','log'])
    except:
        print(f'Usage: {sys.argv[0]} directory [csv|log]', file=sys.stderr)
        sys.exit(1)

    cols = Columns(R,L)
    cols.hline()
    cols('iteration count', 'file')
    tot = 0
    for pwd, _, fns in os.walk(dir_):
        cols.hline()
        cols('', pwd + '/')
        cols.hline('-')

        for fn in sorted(fns):
            if fn.endswith('.log') and ext == 'log':
                with open(os.path.join(pwd, fn), 'r') as f:
                    cur = 0
                    for line in f:
                        try:
                            fields = line.split()
                            if fields[0] == 'its:':
                                cur += int(fields[1])
                        except:
                            pass
                    cols(nice_int(cur) if cur > 0 else '-', fn)
                    tot += cur
            if fn.endswith('.csv') and ext == 'csv':
                with open(os.path.join(pwd, fn), 'r') as f:
                    cur = 0
                    for line in f:
                        try:
                            fields = line.split(',')
                            cur += int(fields[2])
                        except:
                            pass
                    cols(nice_int(cur) if cur > 0 else '-', fn)
                    tot += cur

    cols.hline()
    cols(nice_int(tot), 'TOTAL')
    cols.hline()
    print(cols)
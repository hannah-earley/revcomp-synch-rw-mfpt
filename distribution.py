#!/usr/bin/env python3
import argparse
import sys
import common

class Histogram:
    def __init__(self, binning=1):
        self.binning = binning
        self.counts = {}
        self.total = 0

    def record_event(self, n, count=1):
        self[n] += count
        self.total += count

    def which_bin(self, n):
        binn = self.binning
        return (n // binn) * binn

    def load(self, file):
        self._load(common.LineIterator(file))

    def __iter__(self):
        binn = self.binning
        cl = list(self.counts.items())
        cl.sort()
        for r,c in cl:
            yield r, r + binn, c

    def columns(self):
        return zip(*self)

    def _load(self, lines):
        for line in lines:
            row0, rowf, count = [int(x.strip()) for x in line.split(',')]

            if self.binning is None:
                assert not self.counts and not self.total
                # guess binning...
                binn_ = rowf - row0 + 1
                offs = row0 % binn_
                if offs != 0:
                    raise ValueError("Input file has invalid binning! Found %d (offset %d)..." % (binn_, offs))
                self.binning = binn_
            
            b0 = self.which_bin(row0)
            bf = self.which_bin(rowf)
            if b0 != bf:
                binn = self.binning
                binn_ = rowf - row0 + 1
                offs = row0 % binn
                raise ValueError("Input file has incompatible binning! Found %d (offset %d), expecting %d (offset 0)..." % (binn_, offs, binn))

            self.record_event(row0, count)

    def __getitem__(self, n):
        b = self.which_bin(n)
        return self.counts.setdefault(b, 0)

    def __setitem__(self, n, c):
        b = self.which_bin(n)
        self.counts[b] = c

    def dump_csv(self, file=sys.stdout):
        print("#row0,rowf,count", file=file)
        for r0,rf,c in self:
            print("%d,%d,%d" % (r0,rf-1,c), file=file)

class ExactHist:
    def raw(self, row):
        return 0
    def between(self, row0, rowf):
        return sum(self.raw(row) for row in range(row0, rowf))

class Exact1D:
    class Standard(ExactHist):
        def __init__(self, bias):
            self.bias = bias
            self.p = 0.5 * (1 + bias)
            self.q = 1 - self.p
            self.x0 = bias / self.p
            self.t = (1 - bias) / (1 + bias)

        def raw(self, row):
            row = abs(row)
            return self.x0 * (self.t ** row)

    class Teleporting(Standard):
        def __init__(self, bias, distance):
            super().__init__(bias)
            if bias == 0:
                return

            self.distance = s = distance
            del self.x0
            self.pb = self.p / bias
            self.x1 = 1 / (self.pb * s)
            self.x1pb = 1 / s
            self.xs = self.x1pb * (1 - self.t**s)

        def raw(self, row):
            if self.bias == 0:
                return 0
            
            row = abs(row)
            if row >= self.distance:
                return self.xs * (self.t ** (row - self.distance))
            return self.x1pb * (1 - self.t**row)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-b', '--binning', default=1, type=int)
    parser.add_argument('-s', '--skip', default=0, type=int)
    parser.add_argument('-n', '--limit', default=float('inf'), type=int)
    parser.add_argument('-v', '--every', default=None, type=int)
    parser.add_argument('file', nargs='?')
    args = parser.parse_args()

    skip = args.skip
    binn = args.binning
    file = args.file
    every = args.every
    limit = args.limit

    hist = Histogram(binn)

    if file:
        try:
            with open(file, 'r') as fh:
                hist.load(fh)
        except FileNotFoundError:
            pass

    if limit:
        try:
            lines = common.SkipIterator(sys.stdin, skip)
            for line in common.ProgressIterator(lines, limit, every):
                pos, *_ = [int(x.strip()) for x in line.split(',')]
                hist.record_event(-pos)
        except KeyboardInterrupt:
            pass

    if file:
        with open(file, 'w') as fh:
            hist.dump_csv(fh)
    else:
        hist.dump_csv()
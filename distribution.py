#!/usr/bin/env python3
import argparse
import sys
import common

class Histogram:
    def __init__(self, binning=1):
        self.binning = binning
        self.counts = {}
        self.total = 0
        pass

    def record_event(self, n, count=1):
        self[n] += count
        self.total += count

    def which_bin(self, n):
        binn = self.binning
        return (n // binn) * binn

    def load(self, file):
        self._load(common.LineIterator(file))

    def _load(self, lines):
        for line in lines:
            row0, rowf, count = [int(x.strip()) for x in line.split(',')]
            
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
        cl = list(self.counts.items())
        cl.sort()
        binn = self.binning

        print("#row0,rowf,count", file=file)
        for r,c in cl:
            print("%d,%d,%d" % (r,r+binn-1,c), file=file)


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
#!/usr/bin/env python3
import argparse
import sys

def output_csv(counts, file=sys.stdout):
    cl = list(counts.items())
    cl.sort()

    print("#row0,rowf,count", file=file)
    for r,c in cl:
        print("%d,%d,%d" % (r,r+binn-1,c), file=file)

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-b', '--binning', default=1, type=int)
    parser.add_argument('-B', action='store_true',
                        help='rebin input file if possible')
    parser.add_argument('-s', '--skip', default=0, type=int)
    parser.add_argument('-n', '--limit', default=float('inf'), type=int)
    parser.add_argument('-v', '--every', default=float('inf'), type=int)
    parser.add_argument('-V', action='store_true', help='every 5%%')
    parser.add_argument('file', nargs='?')
    args = parser.parse_args()

    skip = args.skip
    binn = args.binning
    rebin = args.B
    file = args.file
    every = args.every
    limit = args.limit
    n = limit + 0
    m = 0
    verb_ = False

    if args.V:
        every = limit / 20

    counts = {}

    if file:
        try:
            with open(file, 'r') as fh:
                for line in fh:
                    line = line.strip()
                    if line[0] == '#':
                        continue
                    row0, rowf, count = [int(x.strip()) for x in line.split(',')]
                    binn_ = rowf - row0 + 1
                    if binn_ != binn:
                        if rebin and binn % binn_ == 0:
                            row0 = (row0 // binn) * binn
                        else:
                            raise ValueError("Input file has incorrect binning! Found %d, expecting %d..." % (binn_, binn))
                    counts.setdefault(row0, 0)
                    counts[row0] += count
        except FileNotFoundError:
            pass

    if n:
        try:
            for line in sys.stdin:
                if skip > 0:
                    skip -= 1
                    continue

                n -= 1
                m += 1
                pos, *_ = [int(x.strip()) for x in line.split(',')]
                row = ((-pos) // binn) * binn

                counts.setdefault(row, 0)
                counts[row] += 1

                if m % every == 0:
                    if limit < float('inf'):
                        stat = '%.1f%%..' % (100*m/limit)
                    else:
                        stat = '%d..' % m
                    print(stat, file=sys.stderr, end='', flush=True)
                    verb_ = True
                if n <= 0:
                    break
        except KeyboardInterrupt:
            pass

    if verb_:
        print('', file=sys.stderr)

    if file:
        with open(file, 'w') as fh:
            output_csv(counts, fh)
    else:
        output_csv(counts)
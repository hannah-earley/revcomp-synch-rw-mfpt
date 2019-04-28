#!/usr/bin/env python3
import os.path
import argparse
import re

open_mode = 'x'
mean_re = re.compile(r'^mean:([0-9.e+-]+) \(±([0-9.e+-]+)\) var:([0-9.e+-]+)$')

def dat_count(dat, n_=None):
    n = 0
    eof = '# eof'

    try:
        with open(dat, 'r') as f:
            for line in f:
                line = line.strip()
                if line[0] == '#':
                    if line == eof:
                        return n
                    continue
                if line:
                    n += 1
    except FileNotFoundError:
        print('Data file not found: %s' % dat)
        return n_

    if n != n_:
        print("Possibly corrupt: %s (%d)" % (dat,n))
    return n

def rawfine(fin, fout, dat=None):
    n = None
    m = None
    s = None

    if dat is not None:
        n = dat_count(dat, n)

    try:
        with open(fin, 'r') as hin, open(fout, open_mode) as hout:
            hout.write('#mean,stderr,weight\n')

            for line in hin:
                line = line.strip()

                fields = [x.strip() for x in line.split(':')]
                if len(fields) == 2:
                    key, val = fields
                    if key == 'Walker Count':
                        n = int(val)
                    if key == 'Measurement Count':
                        m = int(val)
                    if key == 'Sample Window':
                        s = int(val)
                    if key == 'Persisting to':
                        dat = val
                        n = dat_count(dat, n)

                if line.startswith('mean:'):
                    matches = mean_re.match(line)
                    its = n * m * s
                    mean = float(matches.group(1))
                    err = float(matches.group(2))
                    var = float(matches.group(3))
                    mean_ = 1 / mean
                    err_ = err * mean_ * mean_
                    hout.write('%r,%r,%d\n' % (mean_, err_, its))

    except FileExistsError:
        print('%s: Refusing to overwrite %s' % (fin, fout))

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="")
    parser.add_argument('-x', default='csv~',
                        help='new extension')
    parser.add_argument('-f', default=False, action='store_true',
                        help='overwrite if exists')
    parser.add_argument('file', nargs='*',
                        help='(dirty) log file')

    args = parser.parse_args()
    ext = args.x

    if args.f:
        open_mode = 'w'

    for fin in args.file:
        fbase,_ = os.path.splitext(fin)
        fout = fbase + '.' + ext
        fdat = fbase + '.dat'
        rawfine(fin, fout, fdat)
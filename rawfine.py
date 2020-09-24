#!/usr/bin/env python3
import os.path
import argparse
import datetime
import re
import math

open_mode = 'x'
its_re = re.compile(r'^its: +([0-9]+)$')
mean_re = re.compile(r'^mean:([0-9.e+-]+|inf) \(±([0-9.e+-]+|nan)\) var:([0-9.e+-]+|nan)$')
syst_re = re.compile(r'^(CPU|Wall): +(.*) \((.*)\)$')
real_re = re.compile(r'^Real: (.*) - (.*)$')
real_fmt = '%H:%M:%S %d/%m/%y'
unit_durs = {'y':365*86400,'d':86400,'h':3600,'m':60,'s':1,
             'ms':1e-3,'us':1e-6,'ns':1e-9,'ps':1e-12,
             'fs':1e-15,'as':1e-18,'zs':1e-21,'ys':1e-24}
sims = {"MFPT - 1D Walk": '1d'
       ,"MFPT - 2D Walk (Constrained/Quadrant)": '2d'
       ,"Unit Tests": 'tu'
       ,"Test Bed...": 'tb'}

DAT_CACHE = {}
DAT_EOF = '# eof'
def dat_count(dat, n_=None):
    poss_corrupt = True

    if dat in DAT_CACHE:
        n, poss_corrupt = DAT_CACHE[dat]

    else:
        n = 0
        try:
            with open(dat, 'r') as f:
                for line in f:
                    line = line.strip()
                    if line[0] == '#':
                        if line == DAT_EOF:
                            poss_corrupt = False
                            break
                        continue
                    if line:
                        n += 1
        except FileNotFoundError:
            print('Data file not found: %s' % dat)
            n = None
        DAT_CACHE[dat] = n, poss_corrupt

    if n is None:
        print('Data file not found: %s' % dat)
        DAT_CACHE.pop(dat, None)
        return n_
    if n != n_ and poss_corrupt:
        print("Possibly corrupt: %s (%d,%r)" % (dat,n,n_))
        raise ValueError
    return n

def parse_duration(s):
    toks = s.split()
    vals = toks[0::2]
    units = toks[1::2]
    return sum(float(v) * unit_durs[u] for v,u in zip(vals, units))

def read_log(fin, dat=None, meta=True):
    outp = {
        'sim':None,
        'bias':None, 'dist':None, 'width':None,
        'n':None, 'm':None, 's':None, 'its':None,
        'mean':None, 'var':None, 'err':None,
        'dur':None, 'prog':None,
        'wall':(None,None), 'cpu':(None,None)
    }
    outq = None

    if dat is not None:
        outp['n'] = dat_count(dat, outp['n'])

    try:
        with open(fin, 'r') as hin:

            for line in hin:
                line = line.strip()

                if line.count(':') == 1:
                    key, val = line.split(':', 1)
                    key = key.strip()
                    val = val.strip()

                    if key == 'Running simulation':
                        outp['sim'] = sims[val]
                    if key == 'Bias':
                        outp['bias'] = float(val)
                    if key == 'Distance':
                        outp['dist'] = int(val)
                    if key == 'Constriction Width':
                        outp['width'] = int(val)
                    if key == 'Walker Count':
                        outp['n'] = int(val)
                    if key == 'Measurement Count':
                        outp['m'] = int(val)
                    if key == 'Sample Window':
                        outp['s'] = int(val)
                    if key == 'Persisting to':
                        if outp['n'] is None:
                            dat = val
                            outp['n'] = dat_count(dat, outp['n'])

                if line.startswith('its:'):
                    matches = its_re.match(line)
                    outp['its'] = int(matches.group(1))

                if line.startswith('mean:'):
                    matches = mean_re.match(line)
                    if not matches:
                        print(fin, line)
                        raise ValueError

                    mean_ = float(matches.group(1))
                    err_ = float(matches.group(2))
                    var_ = float(matches.group(3))

                    fin = math.isfinite
                    if not (fin(mean_) and fin(err_) and fin(var_)):
                        if outq:
                            yield outq
                        outq = None
                        continue

                    outp['mean'] = mean_
                    outp['err'] = err_
                    outp['var'] = var_

                    if outq:
                        yield outq
                    outq = outp.copy()
                    outp['its'] = None

                elif line.startswith('Real:'):
                    if meta:
                        matches = real_re.match(line)
                        t0 = datetime.datetime.strptime(matches.group(1), real_fmt)
                        t1 = datetime.datetime.strptime(matches.group(2), real_fmt)
                        outp['dur'] = t1 - t0
                        if outq is not None:
                            outq['dur'] = outp['dur']

                elif line.startswith('CPU:') or line.startswith('Wall:'):
                    if meta:
                        matches = syst_re.match(line)
                        type_ = matches.group(1).lower()
                        type__ = type_ + '_'
                        loop_ = matches.group(2)
                        tot_ = matches.group(3)
                        loop = parse_duration(loop_)
                        tot = parse_duration(tot_)

                        outp[type_] = (loop, tot)
                        outp[type__] = (loop_, tot_)
                        if outq is not None:
                            outq[type_] = outp[type_]
                            outq[type__] = outp[type__]

                elif line.startswith('0%..'):
                    progs = line.split('..')
                    try:
                        if   progs[-1]: outp['prog'] = progs[-1]
                        elif progs[-2]: outp['prog'] = progs[-2]
                    except IndexError:
                        pass

            if outq:
                yield outq
            if outp != outq and outq is not None:
                o1 = outp.copy()
                o2 = outq.copy()
                o1['its'] = o2['its']
                o1['prog'] = o2['prog']
                if o1 != o2:
                    yield outp

    except FileExistsError:
        print('%s: Refusing to overwrite %s' % (fin, fout))

def log2csv(outps):
    for result in outps:
        mean = 1 / result['mean']
        err = result['err'] * mean * mean
        its = result['its']
        if its is None:
            its = result['n'] * result['m'] * result['s']
        yield mean, err, its

def read2csv(*args, **kwargs):
    yield from log2csv(read_log(*args, **kwargs))

def rawfine(fin, fout, dat=None):
    try:
        with open(fout, open_mode) as hout:
            hout.write('#mean,stderr,weight\n')
            for mean, err, its in read2csv(fin, dat):
                hout.write('%r,%r,%d\n' % (mean, err, its))

    except FileExistsError:
        print('%s: Refusing to overwrite %s' % (fin, fout))

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="")
    parser.add_argument('-x', default='csv~',
                        help='new extension')
    parser.add_argument('-f', default=False, action='store_true',
                        help='overwrite if exists')
    parser.add_argument('--dry', default=False, action='store_true')
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

        if args.dry:
            for result in read_log(fin, fdat):
                print(result)
        else:
            rawfine(fin, fout, fdat)
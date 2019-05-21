#!/usr/bin/env python3
import os
import json
import pathlib
import argparse
import rawfine

import common
from common import cmp, py3_cmp, cmp_float

COMMENT_LINE = '#'

EXTS = {
    'log': '.log',
    'persist': '.dat',
    'persist-bak': '.dat~',
    'outp': '.csv',
    'outp_raw': '.csv~',
    'outp_fix': '.csv~~',
    'distr': '.dist',
    'override': '.override'
}

EXTS_rev = {v:k for k,v in EXTS.items()}

TOLERANCE = {
    'bias': 1e-5,
    'distance': 1e-3,
    'datum': 1e-10
}

def handler_common(args):
    try:
        if args.index:
            handler_index(args)
    except AttributeError:
        pass

walk_args = argparse.ArgumentParser()
walk_args.add_argument('-1', dest='sim', action='store_const', const='1d')
walk_args.add_argument('-2', dest='sim', action='store_const', const='2d')
walk_args.add_argument('-b', dest='bias', type=float)
walk_args.add_argument('-d', dest='dist', type=int)
walk_args.add_argument('-w', dest='width', type=int)
walk_args.add_argument('-v', dest='verbose', action='store_true')
def walk_parse(args):
    if isinstance(args, str):
        args = args.split()
    known, unknown = walk_args.parse_known_args(args)
    return known

class Parameters:
    def __init__(self, sim, bias=0, dist=1, width=1):
        self.simulation = sim
        self.bias = float(bias)
        self.distance = int(dist)
        self.width = None
        if width is not None:
            self.width = int(width)

    def __repr__(self):
        return repr((self.simulation, self.bias, self.distance, self.width))

    def __eq__(self, other):
        if self.simulation != other.simulation:
            return False

        if cmp_float(self.bias, other.bias, TOLERANCE['bias'], True) != 0:
            return False
        if cmp_float(self.distance, other.distance, TOLERANCE['distance'], True) != 0:
            return False

        if self.simulation == '2d' and self.width != other.width:
            return False

        return True

class Datum(py3_cmp):
    def __init__(self, mean, stderr, n):
        self.mean = mean
        self.stderr = stderr
        self.n = n

    def __cmp__(self, other):
        cn = cmp(self.n, other.n)
        cm = cmp_float(self.mean, other.mean, TOLERANCE['datum'])
        cs = cmp_float(self.stderr, other.stderr, TOLERANCE['datum'])

        if cn != 0: return cn
        if cm != 0: return cm
        return cs

    @classmethod
    def from_tuple(cls, t):
        return cls(*t)

    @classmethod
    def from_tuples(cls, ts):
        return [cls(*t) for t in ts]

    def __repr__(self):
        return repr((self.mean, self.stderr, self.n))

    def __str__(self):
        return '%r,%r,%d' % (self.mean, self.stderr, self.n)

class Job:
    def __init__(self, jobdir, jobn):
        self.jobdir = jobdir
        self.name = jobn
        self.path = os.path.join(self.jobdir, self.name)
        self.data = set()
        self.unrecognised = []
        self._params = None
        self.overrides = {
            'outp': set()
        }

    def add_data(self, suffix):
        try:
            self.data.add(EXTS_rev[suffix])
        except KeyError:
            self.unrecognised.append(suffix)

    def __repr__(self):
        return 'Job:' + repr((self.jobdir,self.name,self.data,self._params))

    def get_overrides(self):
        if 'override' not in self.data:
            return

        over = self.path + EXTS['override']
        unknown = []
        with open(over, 'r') as f:
            for override in f:
                override = override.strip()
                try:
                    key, arg = override.split(':',1)
                    if key == 'outp':
                        for idx in arg.split(','):
                            if '-' in idx:
                                a,b = map(int, idx.split('-'))
                                for i in range(a,b+1):
                                    self.overrides['outp'].add(int(i))
                            else:
                                self.overrides['outp'].add(int(idx))
                    else:
                        unknown.append(override)
                except ValueError:
                    unknown.append(override)
        if unknown:
            warn = "get_overrides: Job %s: unknown overrides:\n" % self.name
            for over in unknown:
                warn += "  %s\n" % over
            raise Warning(warn)

    def get_parameters(self):
        def method_log():
            dat = None
            if 'persist' in self.data:
                dat = self.path + EXTS['persist']
            if 'log' in self.data:
                fin = self.path + EXTS['log']
                outp = None
                for outp in rawfine.read_log(fin, dat):
                    pass
                if outp is None:
                    return

                sim = outp['sim']
                if sim:
                    return Parameters(sim, outp['bias'],
                                           outp['dist'],
                                           outp['width'])

        def method_dat():
            if 'persist' in self.data:
                dat = self.path + EXTS['persist']
                with open(dat, 'r') as fh:
                    for line in fh:
                        assert line.startswith(COMMENT_LINE)
                        cmd = line[1:].strip()
                        args = walk_parse(cmd)

                        sim = args.sim
                        if sim:
                            return Parameters(sim, args.bias,
                                                   args.dist,
                                                   args.width)
                        break

        def method_name():
            pure_name = pathlib.PurePath(self.name).parts[-1]
            try:
                sim, bias, *rest = pure_name.split('-')
            except ValueError:
                return

            if sim == '1d':
                dist, *_ = rest
                return Parameters(sim, bias, dist)
            elif sim == '2d':
                width, dist, *_ = rest
                return Parameters(sim, bias, dist, width)

        if self._params is None:
            self._params = False
            methods = [method_log, method_dat, method_name]
            for method in methods:
                ret = method()
                if ret:
                    self._params = ret
                    break

        if self._params:
            return self._params
        raise Warning("Couldn't infer job parameters for %s!" % self.name)

    def check_outp(self):
        if 'log' not in self.data:
            raise Warning("check_outp: Job %s missing log file" % self.name)
        if 'outp' not in self.data:
            raise Warning("check_outp: Job %s missing outp file" % self.name)

        log = self.path + EXTS['log']
        outp = self.path + EXTS['outp']
        dat = None

        if 'persist' in self.data:
            dat = self.path + EXTS['persist']

        raw_dat = Datum.from_tuples(rawfine.read2csv(log, dat))
        csv_dat = []

        with open(outp, 'r') as fh:
            for line in fh:
                line = line.strip()
                if not line:
                    continue
                if line.startswith(COMMENT_LINE):
                    continue

                mean, err, its = line.split(',')
                csv_dat.append(Datum(float(mean), float(err), int(its)))

        lcs = common.LCS(csv_dat, raw_dat)
        com_dat = lcs.common()

        diff = lcs.diff()[0]
        if diff:
            fix = "#mean,stderr,weight\n"
            for d,k,x in diff:
                fix += str(x) + "\n"

            warn = "check_outp: Job %s:\n" % self.name
            warn += "  Inconsistent outp(%d, %s) log(%d, %s) outputs!\n" % (
                        len(csv_dat),EXTS['outp'],len(raw_dat),EXTS['log'])

            overs = set()
            notcsv = 0
            warn2 = "    in outp but not log (*):\n"
            for d,k,x in diff:
                if d < 0:
                    if k in self.overrides['outp']:
                        overs.add(k)
                    else:
                        warn2 += "      %4d: %s\n" % (k,x)
                        notcsv += 1
            if notcsv:
                warn2 += "    * may indicate incorrect weighting...\n"
                warn += warn2
            excess = self.overrides['outp'].difference(overs)
            if excess:
                warn += "    overriden unnecessarily: %r\n" % excess

            notlog = 0
            warn2 = "    in log but not outp:\n"
            for d,k,x in diff:
                if d > 0:
                    warn2 += "      %4d: %s (mean~%f)\n" % (k,x,1/x.mean)
                    notlog += 1
            if notlog:
                warn += warn2

            warn += "  N.B. Record indices may not match line numbers\n"
            warn += "  A candidate fix has been written to %s\n" % EXTS['outp_fix'] 

            if notcsv or notlog or excess:
                with open(self.path + EXTS['outp_fix'], 'w') as f:
                    f.write(fix)
                raise Warning(warn)


class Experiment:
    def __init__(self, params):
        self.params = params
        self.jobs = {}

    def add(self, job):
        self.jobs[job.name] = job

    def __repr__(self):
        return repr((self.params, self.jobs))

    def resolve(self):
        raise NotImplementedError()
        # if one job, easy to resolve
        # otherwise, ask what to do...
        pass

    def mfpt(self):
        pass

    def distribution(self):
        pass

    def __iter__(self):
        return iter(self.jobs.items())

class ExperimentEnsemble:
    """
    Experiments are parameterised by:
    - sim (1d or 2d)
    - bias
    - distance
    - [width]

    Cannot use simple dictionary because:
    - bias is floating
    - width is irrelevant for 1d, ignore

    To address, we use a simple association list...
    """

    def __init__(self):
        self.experiments = []

    def find(self, params1):
        for i, (params2, _) in enumerate(self.experiments):
            if params1 == params2:
                return i

        # expt not found, create new expt
        expt = Experiment(params1)
        self.experiments.append((params1, expt))
        return len(self.experiments) - 1

    def __getitem__(self, params):
        index = self.find(params)
        _, expt = self.experiments[index]
        return expt

    def __repr__(self):
        return repr(self.experiments)

    def __iter__(self):
        return iter(self.experiments)

@common.run_once(error=False, warn=False)
def handler_index(args):
    jobdir = args.jobdir
    idxfile = args.index_file
    idxpath = os.path.join(jobdir, idxfile)

    index = {}
    try:
        with open(idxpath, 'r') as idx:
            index = json.load(idx)
    except FileNotFoundError:
        pass

    jobs = {}

    for root,_,files in os.walk(jobdir):
        for file in files:
            path = os.path.join(root, file)
            rpath = pathlib.PurePath(path).relative_to(jobdir)
            jobn = str(rpath.with_suffix(''))

            jobs.setdefault(jobn, Job(jobdir, jobn))
            jobs[jobn].add_data(rpath.suffix)

    expts = ExperimentEnsemble()

    warnings = False
    for jobn,job in jobs.items():
        try:
            job.get_overrides()
            job.check_outp()
            params = job.get_parameters()
        except Warning as w:
            print(w)
            warnings = True
            continue
        expts[params].add(job)

    if warnings:
        print("\n*** Please fix the above warnings before continuing")
        print("* For inconsistent outp/log, manually edit relevant files to resolve inconsistencies or use override*")
        print("* For missing outp files, use rawfine.py cautiously to reconstruct")
        print("* For missing log files, first touch the log file, then override* the relevant warnings.")
        print("* For uninferred parameters, you may have an empty log file - add params to the log file...")
        print("** .override files: see README.md")
        return

    for p,xs in expts:
        print(p,'::')
        for j,x in xs:
            print('  ',j)
            print('   ',x)

    return
    with open(idxpath, 'w') as idx:
        json.dump(index, idx)

def handler_hist(args):
    pass
def handler_plot(args):
    pass

if __name__ == '__main__':
    parser = common.CommandParser(description='Suite of analysis tools for data generated by other tools in this repository.')
    parser.add_argument('-j', '--jobdir', default='./jobs', help='data directory')
    parser.add_argument('-i', '--index', default=False, action='store_true',
        help='automatically reindex first')
    parser.add_argument('--index-file', default='_index.json')



    parser_idx = parser.add_command('index', aliases=['idx'],
                help='Index files for use by other commands.',
                description='Index the files in the data directory by -bdw parameters for use by other analysis tools. For most output files, this process will succeed automatically, but where there are multiple jobs for a given parameter set you will be asked to manually resolve the ambiguity.')
    parser_idx.handler(handler_index)



    parser_hist = parser.add_command('histogram', aliases=['hist'],
                help='Generate histograms / approximate PDF plots.',
                description='Reads in data from distribution.py and produces a plot of the approximate PDF.')
    parser_hist.add_argument('-x', action='store_true',
                help='Superimpose the exact equilibrium distribution for a reflecting boundary condition.')
    # parser_hist.handler(handler_hist)



    parser_plot = parser.add_command('plot',
                help='Generate MFPT plots.',
                description='',
                epilog='RANGE: Can be a comma delimited list, a hyphenated range, or a mix, e.g. "1,3,6"; "5-9"; "1,3-10,4"; etc. If none or "-", then all matching parameters will be selected.')
    parser_plot.add_argument('-T', '--mfpt', action='store_true',
                help='Generate raw mfpt/distance plot')
    parser_plot.add_argument('-R', '--residue', action='store_true',
                help='Generate residue plot to estimate penalty')
    parser_plot.add_argument('-u', type=float, default=1.5,
                help='Coefficient for unpenalised MFPT used to compute residues.')
    parser_plot.add_argument('-s', '--series', choices="bdw",
                help='Superimpose the specified data series on each plot.')
    parser_plot.add_argument('-b', '--biases', metavar='RANGE', nargs='?',
                help='Which bias(es) to plot.')
    parser_plot.add_argument('-d', '--distances', metavar='RANGE', nargs='?',
                help='Which distance(s) to plot.')
    parser_plot.add_argument('-w', '--widths', metavar='RANGE', nargs='?',
                help='Which width(s) to plot.')
    parser_plot.add_argument('-f', '--format', choices=['pgf','png','eps','gui'], default='gui')
    # parser_plot.handler(handler_plot)



    args = parser.parse_args()
    handler_common(args)
    args.handler(args)

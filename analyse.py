#!/usr/bin/env python3
import os
import json
import pathlib
import argparse
from functools import wraps
import itertools as it
import stat
import math

import refine
import rawfine
import common
import distribution
from common import cmp, py3_cmp, cmp_float
common.SHOW_TIMERS = False

import matplotlib
import matplotlib.pyplot as plt
import numpy as np

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

EXT_PLAN = '.plan'

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
walk_args.add_argument('-b', dest='bias', type=float, default=0)
walk_args.add_argument('-d', dest='dist', type=int, default=1)
walk_args.add_argument('-w', dest='width', type=int, default=1)
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

    def __str__(self):
        sim = self.simulation
        bias = self.bias
        dist = self.distance
        width = self.width

        if sim == '1d':
            return '%s-%r-%d' % (sim, bias, dist)
        elif sim == '2d':
            return '%s-%r-%d-%d' % (sim, bias, width, dist)

        raise Warning('Parameters.__str__: cannot convert simulation type "%s" to a path name...' % sim)

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

    @staticmethod
    def from_log(fin, dat=None):
        outp = None
        for outp in rawfine.read_log(fin, dat, meta=False):
            pass
        if outp is None:
            return

        sim = outp['sim']
        if sim:
            return Parameters(sim, outp['bias'],
                                   outp['dist'],
                                   outp['width'])

    @staticmethod
    def from_dat(dat):
        with open(dat, 'r') as fh:
            for line in fh:
                assert line.startswith(COMMENT_LINE)
                cmd = line[1:].strip()
                return Parameters.from_args(cmd)

    @staticmethod
    def from_args(cmd):
        args = walk_parse(cmd)
        sim = args.sim
        if sim:
            return Parameters(sim, args.bias,
                                   args.dist,
                                   args.width)


    @staticmethod
    def from_name(pure_name):
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

    def path_for(self, type_, force=False):
        if force or type_ in self.data:
            return self.path + EXTS[type_]

    def add_data(self, suffix):
        try:
            self.data.add(EXTS_rev[suffix])
        except KeyError:
            self.unrecognised.append(suffix)

    def __repr__(self):
        return 'Job:' + repr((self.jobdir,self.name,self.data,self._params))

    def get_overrides(self):
        over = self.path_for('override')
        if not over:
            return

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
            dat = self.path_for('persist')
            fin = self.path_for('log')
            if fin:
                return Parameters.from_log(fin, dat)

        def method_dat():
            dat = self.path_for('persist')
            if dat:
                return Parameters.from_dat(dat)

        def method_name():
            pure_name = pathlib.PurePath(self.name).parts[-1]
            return Parameters.from_name(pure_name)

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
        log = self.path_for('log')
        outp = self.path_for('outp')
        dat = self.path_for('persist')

        if not log:
            raise Warning("check_outp: Job %s missing log file" % self.name)
        if not outp:
            raise Warning("check_outp: Job %s missing outp file" % self.name)

        with common.timer('rawfine'):
            raw_dat = Datum.from_tuples(rawfine.read2csv(log, dat, meta=False))

        with common.timer('refine'):
            csv_dat = []
            with open(outp, 'r') as fh:
                for line in fh:
                    line = line.strip()
                    if not line:
                        continue
                    if line.startswith(COMMENT_LINE):
                        continue

                    mean, err, its = line.split(',')
                    mean = float(mean)
                    err = float(err)
                    its = int(its)
                    fin = math.isfinite
                    if not (fin(mean) and fin(err)):
                        continue

                    csv_dat.append(Datum(mean, err, its))

        # lcs = common.LCS(csv_dat, raw_dat)
        with common.timer('lcs'):
            lcs = common.LCSGreedyApprox(raw_dat, csv_dat)
        lcs.invert()
        com_dat = lcs.common()

        diff = lcs.diff()[0]
        if diff:
            fix = "#mean,stderr,weight\n"
            for d,k,x in diff:
                try:
                    x,_ = x
                except TypeError:
                    pass
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
                with open(self.path_for('outp_fix', True), 'w') as f:
                    f.write(fix)
                raise Warning(warn)


class Experiment:
    def __init__(self, args, params):
        self.args = args
        self.params = params
        self.jobs = {}
        self.plan = None

    def add(self, job):
        self.jobs[job.name] = job

    def __repr__(self):
        return repr((self.params, self.jobs))

    name_cache = None
    @classmethod
    def scan_names(cls, jobdir):
        if cls.name_cache is not None:
            return
        cls.name_cache = []

        for root,_,files in os.walk(jobdir):
            for file in files:
                path = os.path.join(root, file)
                ppath = pathlib.PurePath(path)
                rpath = ppath.relative_to(jobdir)
                if rpath.suffix != EXT_PLAN:
                    continue

                name = rpath.with_suffix('')
                bpath = ppath.with_suffix('')
                pname = name.name
                params = Parameters.from_name(name.name)

                if params:
                    for name2,_,params2 in cls.name_cache:
                        if params == params2:
                            warn = "Experiment.scan_names: duplicate plans found:\n"
                            warn += "  %s :: %s\n" % (name2, params2)
                            warn += "  %s :: %s\n" % (name, params)
                            cls.name_cache = None
                            raise Warning(warn)
                    cls.name_cache.append((str(name), str(bpath), params))

    def name(self):
        """Generate a name for this experiment, which may consist of a group of jobs (and may acquire more jobs over time. This presents a couple problems:
        - the bias parameter is a float, and so doesn't necessarily have a best unique representation
        - we have some tolerance on bias and distance, so adding more jobs over time could perturb which parameters we use to build the name
        - in addition, the order of jobs is not necessarily stable so even when not adding jobs the parameters could change order...

        To address these, we should first search the directory for names whose parameters are compatible with ours. If none exists, create a new one from our parameter set. If the user doesn't like the chosen name, they can perturb it to a similar but consistent name.

        We should also cache the names so each experiment doesn't rescan the directory"""

        self.scan_names(self.args.jobdir)
        params = self.params
        for name2,path2,params2 in self.name_cache:
            if params == params2:
                return name2, path2

        name = str(params)
        path = os.path.join(self.args.jobdir, name)
        self.name_cache.append((name, path, params))
        return name, path

    def path_for(self, ext=EXT_PLAN):
        _, path = self.name()
        return path + ext

    def resolve(self):
        if self.plan is None:
            name, _ = self.name()
            pplan = self.path_for()

            jobs = list(self.jobs.keys())
            jobs.sort()
            plan = {}
            changed = False

            try:
                if self.args.replan:
                    raise FileNotFoundError

                # load plan from last index
                with open(pplan, 'r') as f:
                    plan = json.load(f)

                # update plan if jobset changed
                jobs1 = set(jobs)
                jobs2 = set(j['job'] for j in plan['sequence'])

                for j in plan['sequence']:
                    if j['job'] not in jobs1:
                        j['deleted'] = True
                        plan['staged'] = True
                        changed = True

                for job in jobs:
                    if job not in jobs2:
                        plan['sequence'].append({
                            "job": job,
                            "range": "%d:" % self.args.skip
                        })
                        plan['staged'] = True
                        changed = True

            except (json.decoder.JSONDecodeError, FileNotFoundError, KeyError) as e:
                if not isinstance(e, FileNotFoundError):
                    print("Error encountered when loading experiment plan %s:" % name)
                    print("   " + str(e))
                    print("   ... replanning ...")

                # construct initial plan
                plan = {
                    "staged": True,
                    "sequence": [{
                        "job": job,
                        "range" : "%d:" % self.args.skip
                    } for job in jobs]
                }

                if not self.args.stage_single and len(jobs) == 1:
                    del plan['staged']

                changed = True

            if changed:
                with open(pplan, 'w') as f:
                    json.dump(plan, f, sort_keys=True, indent=4)

            self.plan = plan

        if self.plan.get('staged', False):
            warn = "Experiment.resolve: %s:\n" % self.params
            warn += "  experiment plan currently staged...\n"

            raise Warning(warn)

    @staticmethod
    def read_range(r, type_=lambda x:slice(*x)):
        fields = []
        for field in r.split(':'):
            field = field.strip()
            try:
                field = int(field)
            except ValueError:
                field = None
            fields.append(field)
        return type_(fields)

    def mfpt(self):
        data = refine.Experiment()
        for dataset in self.plan['sequence']:
            range_ = self.read_range(dataset['range'], list)
            job = self.jobs[dataset['job']]
            outp = job.path_for('outp')
            if outp:
                with open(outp, 'r') as f:
                    lines = common.LineIterator(f, COMMENT_LINE)
                    data._update(it.islice(lines, *range_))

        try:
            return data.summarise()
        except ValueError:
            pass

    def mtimes(self):
        def gmt(path, old=None):
            try:
                new = os.path.getmtime(path)
            except (TypeError, FileNotFoundError):
                return old
            if old:
                return max(new, old)
            return new

        tdistr = None
        thist = None

        for dataset in self.plan['sequence']:
            job = self.jobs[dataset['job']]

            distr = job.path_for('distr')
            tdistr = gmt(distr, tdistr)

        hist = self.path_for('.png')
        thist = gmt(hist, thist)

        return {
            'distr': tdistr,
            'hist': thist
        }




    def distribution(self):
        hist = distribution.Histogram(None)
        read = set()

        for dataset in self.plan['sequence']:
            job = self.jobs[dataset['job']]
            distr = job.path_for('distr')
            if distr:
                with open(distr, 'r') as f:
                    stats = os.fstat(f.fileno())
                    fid = stats.st_dev, stats.st_ino
                    if fid not in read:
                        read.add(fid)
                        hist.load(f)

        return hist


    def __iter__(self):
        return iter(self.jobs.items())

class Index:
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

    def __init__(self, args):
        self.args = args
        self.experiments = []

    def find(self, params1, create=True):
        for i, (params2, _) in enumerate(self.experiments):
            if params1 == params2:
                return i

        if create:
            # expt not found, create new expt
            expt = Experiment(self.args, params1)
            self.experiments.append((params1, expt))
            return len(self.experiments) - 1

    def __getitem__(self, params):
        index = self.find(params)
        _, expt = self.experiments[index]
        return expt

    def __contains__(self, params):
        return self.find(params, False) is not None

    def __repr__(self):
        return repr(self.experiments)

    def __iter__(self):
        return iter(self.experiments)

@common.run_once(error=False, warn=False)
def generate_index(args):
    jobdir = args.jobdir
    jobs = {}

    with common.timer(1):
        for root,_,files in os.walk(jobdir):
            for file in files:
                path = os.path.join(root, file)
                rpath = pathlib.PurePath(path).relative_to(jobdir)
                jobn = str(rpath.with_suffix(''))

                if rpath.suffix in EXTS_rev:
                    jobs.setdefault(jobn, Job(jobdir, jobn))
                    jobs[jobn].add_data(rpath.suffix)

    index = Index(args)

    with common.timer(2):
        warnings = False
        for jobn,job in jobs.items():
            print('.', end='', flush=True)
            with common.timer((3,jobn)):
                try:
                    with common.timer((3,jobn,'overrides')):
                        job.get_overrides()
                    with common.timer((3,jobn,'outp')):
                        job.check_outp()
                    with common.timer((3,jobn,'params')):
                        params = job.get_parameters()
                except Warning as w:
                    print(w)
                    warnings = True
                    continue
                index[params].add(job)
        print('')

        if warnings:
            print("\n*** Please fix the above warnings before continuing")
            print("* For inconsistent outp/log, manually edit relevant files to resolve inconsistencies or use override*")
            print("* For missing outp files, use rawfine.py cautiously to reconstruct")
            print("* For missing log files, first touch the log file, then override* the relevant warnings.")
            print("* For uninferred parameters, you may have an empty log file - add params to the log file...")
            print("** .override files: see README.md")
            return False

    with common.timer(4):
        warnings = False
        for p,xs in index:
            try:
                xs.resolve()
            except Warning as w:
                print(w)
                warnings = True
                continue

        if warnings:
            print("\n*** Please fix the above warnings before continuing")
            print("* Staged plans need manual review before use;")
            print("  if the plan looks right, remove the staged key")
            return False

    return index

def indexed_handler(f):
    @wraps(f)
    def wrapped(args):
        index = generate_index(args)
        if index:
            return f(args, index)
        return False
    return wrapped

@indexed_handler
def handler_index(args, index):
    for params, expt in index:
        print(params, '::')
        for jobn,job in expt:
            print('  ',jobn)
            print('   ',job)
        print(' ',expt.plan)

@indexed_handler
def handler_mfpt(args, index):
    for params, expt in index:
        print(params, '::', expt.mfpt())

@indexed_handler
def handler_hist(args, index):
    for params, expt in index:
        mtimes = expt.mtimes()
        if mtimes['distr'] and mtimes['hist']:
            if mtimes['hist'] > mtimes['distr']:
                if not args.all:
                    continue

        hist = expt.distribution()
        file = expt.path_for('.png')

        if hist.total:
            row0s, rowfs, counts = hist.columns()
            freqs = [c/hist.total for c in counts]

            fig, ax = plt.subplots()
            ax.plot(row0s, freqs)

            if params.simulation == '1d':
                ex = distribution.Exact1D.Teleporting(params.bias, params.distance)
                freqs2 = [ex.between(r0,rf) for r0,rf in zip(row0s, rowfs)]
                ax.plot(row0s, freqs2)

            elif params.simulation == '2d':
                ex = distribution.Exact2D.Reversible(params.bias, params.width)
                freqs2 = [ex.between(r0,rf) for r0,rf in zip(row0s, rowfs)]
                ax.plot(row0s, freqs2)

            ax.grid()
            fig.savefig(file)

            if params.simulation == '2d':
                # plt.show()
                pass

            plt.close(fig)

            if params.simulation == '2d':
                file2 = expt.path_for('-col-avg.png')
                widths = (np.array(row0s) + np.array(rowfs)) * 0.5
                afreqs = np.array(freqs) / widths
                fig, ax = plt.subplots()
                ax.plot(row0s, afreqs)
                fig.savefig(file2)
                plt.close(fig)

            print(file)


@indexed_handler
def handler_dist(args, index):
    """
for width in 1 10 100 1000; for dist in 1 3 10 31 100 316 1000 3162 10000 31622 100000 316227 1000000; ./job.sh -b 0.01 -w $width -d $dist -S width-0.01 -- -n 10000000 -v 100000; end; end

for width in 1 10 100 1000; for dist in 1 3 10 31 100 316 1000 3162 10000 31622 100000 316227 1000000; ./job.sh -b 0.1 -w $width -d $dist -S width-test -- -n 10000000 -v 100000; end; end
    """
    pass


@indexed_handler
def handler_plot(args, index):
    pass

if __name__ == '__main__':
    parser = common.CommandParser(description='Suite of analysis tools for data generated by other tools in this repository.')
    parser.add_argument('-j', '--jobdir', default='./jobs', help='data directory')
    parser.add_argument('--skip', default=2,
        help='default number of outputs to skip on each job (when indexing initially)')
    parser.add_argument('--no-stage-single', dest='stage_single', action='store_false',
        help="don't place an initial hold on experiments consisting of single jobs")
    parser.add_argument('--replan', action='store_true',
        help='scrap old experiment plans and reindex')



    parser_idx = parser.add_command('index', aliases=['idx'],
                help='Just index ',
                # help='Index files for use by other commands.',
                description='Index the files in the data directory by -bdw parameters for use by other analysis tools. For most output files, this process will succeed automatically, but where there are multiple jobs for a given parameter set you will be asked to manually resolve the ambiguity.')
    parser_idx.handler(handler_index)


    parser_mfpt = parser.add_command('mfpt', help='Report MFPTs')
    parser_mfpt.handler(handler_mfpt)


    parser_hist = parser.add_command('histogram', aliases=['hist'],
                help='Generate histograms / approximate PDF plots.',
                description='Reads in data from distribution.py and produces a plot of the approximate PDF. By default only regenerates plots where the distribution is newer than the plot.')
    parser_hist.add_argument('-a', '--all', action='store_true',
                help='Generate histograms for all jobs')
    parser_hist.handler(handler_hist)



    parser_dist = parser.add_command('distribution', aliases=['dist'],
                help='Generate distributions.',
                description='For each job found, generate a distribution. By default we skip jobs which already have an associated .dist file, but this can be overriden with --all.')
    parser_dist.add_argument('-a', '--all', action='store_true',
                help='Generate distributions for all jobs')
    parser_dist.add_argument('--regenerate', action='store_true',
                help='Backup old distributions and regenerate')
    parser_dist.add_argument('-n', '--limit', default=10000000,
                help='Datums to sample, default is 1e7.')
    parser_dist.add_argument('params', nargs=argparse.REMAINDER,
                help='Specify parameters for distribution.py')
    parser_dist.handler(handler_dist)



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
    parser_plot.handler(handler_plot)



    args = parser.parse_args()
    handler_common(args)
    args.handler(args)

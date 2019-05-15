#!/usr/bin/env python3
import os
import common
import json
import pathlib
import argparse
import rawfine

COMMENT_LINE = '#'

EXTS = {
    'log': '.log',
    'persist': '.dat',
    'persist-bak': '.dat~',
    'outp': '.csv',
    'outp_raw': '.csv~',
    'distr': '.dist'
}

EXTS_rev = {v:k for k,v in EXTS.items()}

TOLERANCE = {
    'bias': 1e-5,
    'distance': 1e-3
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

        try:
            if abs(self.bias - other.bias)/self.bias > TOLERANCE['bias']:
                return False
        except ZeroDivisionError:
            if abs(other.bias) > TOLERANCE['bias']:
                return False

        try:
            if abs(self.distance - other.distance)/self.distance > TOLERANCE['distance']:
                return False
        except ZeroDivisionError:
            if abs(other.distance) > TOLERANCE['distance']:
                return False

        if self.simulation == '2d' and self.width != other.width:
            return False

        return True

class Job:
    def __init__(self, jobdir, jobn):
        self.jobdir = jobdir
        self.name = jobn
        self.path = os.path.join(self.jobdir, self.name)
        self.data = set()
        self.unrecognised = []
        self._params = None

    def add_data(self, suffix):
        try:
            self.data.add(EXTS_rev[suffix])
        except KeyError:
            self.unrecognised.append(suffix)

    def __repr__(self):
        return 'Job:' + repr((self.jobdir,self.name,self.data,self._params))

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

        raw_dat = list(rawfine.read2csv(log, dat))
        csv_dat = []

        with open(outp, 'r') as fh:
            for line in fh:
                line = line.strip()
                if not line:
                    continue
                if line.startswith(COMMENT_LINE):
                    continue

                mean, err, its = line.split(',')
                csv_dat.append((float(mean), float(err), int(its)))

        print(len(raw_dat), len(csv_dat))
        # print(list(zip(raw_dat,csv_dat)))


class Experiment:
    def __init__(self, params):
        self.params = params
        self.jobs = {}

    def add(self, job):
        self.jobs[job.name] = job

    def __repr__(self):
        return repr((self.params, self.jobs))

    def resolve(self):
        # if one job, easy to resolve
        # otherwise, ask what to do...
        pass

    def mfpt(self):
        pass

    def distribution(self):
        pass

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

            # print(file, path)

    # print(jobs)
    expts = ExperimentEnsemble()

    for jobn,job in jobs.items():
        try:
            job.check_outp()
            params = job.get_parameters()
        except Warning as w:
            print(w)
            continue
        expts[params].add(job)

    print(expts)

    return
    with open(idxpath, 'w') as idx:
        json.dump(index, idx)

def handler_hist(args):
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
    parser_hist.handler(handler_hist)



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



    args = parser.parse_args()
    handler_common(args)
    args.handler(args)

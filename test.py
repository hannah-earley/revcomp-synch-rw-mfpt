#!/usr/bin/env python3
import numpy as np
import scipy
import scipy.special
import scipy.optimize
import math
import matplotlib.pyplot as plt

fs = ['test4d.log','test4c.log','test4b.log']
fs = ['test3d.log','test3c.log','test3b.log']
fs = ['test3f.log']
r0 = 0
lim = 100
lim_ = 120
bias = 0.01
plot = [140,160,180,200,220,240,260,280,300]
tot = 0

show_mids = False
show_steady_state = False
show_det_bal = False
do_fit = True
show_fit_plots = True

## load data

data = {}
for f in fs:
    with open(f) as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            if line[0] == '#':
                continue

            r,_,_,lat = line.split(' ', 3)
            r = int(r)
            cols = list(map(float, lat[1:-1].split(',')))

            if r in data:
                cols = list(map(sum, zip(cols, data[r])))

            data[r] = cols

data = list(data.items())
data.sort()

# s = 0
# for i,row in data:
#     s += sum(row)
# print('sum',s)
# for i,row in data:
#     for j in range(len(row)):
#         row[j] /= s
# print(data[:5])
# 1/0




if show_mids:

    mids = []
    for r,lat in data:
        l = len(lat)
        if l % 2 == 0:
            mids.append((r, 0.5*(lat[l//2-1]+lat[l//2])))
        else:
            mids.append((r, lat[l//2]))

    for r,mid in mids:
        print('%s\t%r' % (r,mid))



if show_steady_state:

    ddat = {}
    for r, lat in data:
        tot += sum(lat)
    for r, lat in data:
        ddat[r] = [x/tot for x in lat]

    p = 0.5 * (1+bias)
    q = 0.5 * (1-bias)
    p2 = 0.5 * p
    q2 = 0.5 * q

    deltas = {}

    dget = lambda r,i: ddat[r][i]
    dget2 = lambda r,i: ddat[r][i] + ddat[r][i+1]

    if r0+1 in ddat and r0+2 in ddat:
        delta0 = p2*dget2(r0+2,0) - (1-p2)*dget(r0+1,0)
        deltaf = p2*dget2(r0+2,-2) - (1-p2)*dget(r0+1,-1)
        deltam = [p2*dget2(r0+2,i) - dget(r0+1,i) for i in range(1,r0+1)]
        deltas[r0+1] = [delta0] + deltam + [deltaf]

    fig, ax = plt.subplots()
    cs = list(range(0,r0+2))
    rele = [d/x for d,x in zip(deltas[r0+1],ddat[r0+1])]
    ax.plot(cs, rele)
    ax.set_title('r0+1=%d'%(r0+1))
    ax.grid()
    plt.show(fig)
    plt.close(fig)

    for r in range(r0+2,max(ddat)):
        if r-1 not in ddat or r not in ddat or r+1 not in ddat:
            continue

        delta0 = q2*dget(r-1,0) + p2*dget2(r+1,0) - (1-p2)*dget(r,0)
        deltaf = q2*dget(r-1,-1) + p2*dget2(r+1,-2) - (1-p2)*dget(r,-1)
        deltam = [q2*dget2(r-1,i-1) + p2*dget2(r+1,i) - dget(r,i) for i in range(1,r)]
        deltas[r] = [delta0] + deltam + [deltaf]

    fig, ax = plt.subplots()
    for r,ds in deltas.items():
        cs = list(range(0,r+1))
        try:
            rele = [d/x for d,x in zip(ds,ddat[r])]
            # print(r,rele)
            ax.plot(cs, rele)
        except ZeroDivisionError:
            pass
        if r % 20 == 0:
            ax.set_title("[%d,%d]" % (r-19,r))
            ax.grid()
            plt.show(fig)
            plt.close(fig)
            fig, ax = plt.subplots()
    ax.grid()
    plt.show(fig)
    plt.close(fig)
    fig, ax = plt.subplots()

if show_det_bal:
    pass




if do_fit:

    ## fit to ansatz (cosh,cos)

    def ansatz1(x, const, mult, weight):
        a = 1
        f = np.cosh
        if weight < 0:
            wight = -weight
            a = -1
            f = np.cos
        return const + mult*a*(f(weight*0.5*x) - 0*1)

    def ansatz2(x, const, mult, weight):
        return const + mult*0.5*(1-np.cos(weight*0.5*x))

    def ansatz3(x, const, mult1, weight1, mult2, weight2):
        return const + mult1*(np.cosh(weight1*0.5*x)-1) + mult2*0.5*(1-np.cos(weight*0.5*x))

    def ansatz4(x, const, mult, weight, a):
        return const + mult*(a*(np.cosh(weight*0.5*x)-1) + (1-a)*(1-np.cos(weight*0.5*x)))

    def ansatz5(r):
        def my_bin(r,c):
            gamma = scipy.special.gamma
            a = gamma(0.5*r)
            b = a/gamma(c)
            c = a/gamma(r-c)
            return b * c

        def ansatz_(x, alph,beta,gamm):
            c = (x+r)//2
            return alph + beta * (1 - my_bin(r/gamm, c/gamm))
        return ansatz_

    def ansatz6(r):
        def my_harm(n):
            return scipy.special.digamma(n+1) + np.euler_gamma

        def my_betabin1(n,k):
            h = my_harm
            hn2 = h(0.5 * n)
            return (h(k)-2*hn2+h(n-k)) / (-2*hn2+h(n))

        def my_betabin(n,k,a):
            if a == 1:
                return my_betabin1(n,k)

            beta = scipy.special.beta
            bino = scipy.special.binom
            n2 = 0.5 * n

            bana = beta(a+n,a)
            deno = bana - beta(a+n2, a+n2)*bino(n, n2)
            numo = bana - beta(a+k, a+n-k)*bino(n, k)
            return 1 - (numo/deno)

        def ansatz_(x, alph,beta,gamm):
            c = (x+r)//2
            return alph + beta*my_betabin(r,c,gamm)

        return ansatz_

    def ansatz7(r):
        def my_betabin(n,k,a):
            beta = scipy.special.beta
            bino = scipy.special.binom
            return bino(n, k)*beta(a+k, a+n-k)/beta(a,a)

        def ansatz_(x, alph,beta,gamm):
            c = (x+r)//2
            return alph - beta*my_betabin(r,c,gamm)

        return ansatz_


    def ansatz1b(w):
        return lambda x,c,m: ansatz1(x,c,m,w)

    rs = []
    consts = []
    mult1s = []
    weight1s = []
    mult2s = []
    weight2s = []
    as_ = []
    errs = []

    for ii,(r,lat) in enumerate(data):
        # if r > lim: break
        # if r < 130: continue
        if r > 600: break

        # ansatz = ansatz1 #3

        xs = np.arange(-r,r+1,2)
        ys = np.array(lat)

        # approximate fit...
        diff = max(ys) - min(ys)
        const = min(ys)
        weight = bias
        # weight = 0.5*np.log((1-bias)/(1+bias))
        # if r > lim:
        #     weight = -weight
        # mult = diff / ansatz1(r,0,1,weight)

        # a = 1 #if r < lim else 0
        # # p0 = (const, mult, weight) #, a) #, (1-a)*mult, bias)
        # p0 = (const, mult)
        # ansatz = ansatz1b(weight)
        # ys0 = ansatz(xs, *p0)

        # ansatz = ansatz5(r)
        # p0 = (const, diff, lim)
        # ys0 = ansatz(xs, *p0)

        # ansatz = ansatz6(r)
        # p0 = (const, diff, 2)
        # ys0 = ansatz(xs, *p0)

        for g0 in [2,0.5]:

            ansatz = ansatz7(r)
            # g0 = 2 #0.5
            min_ = ansatz(0,  0,1,g0)
            max_ = ansatz(r,  0,1,g0)
            mult = diff / (max_-min_)
            min_ = ansatz(0,  0,mult,g0)
            p0 = (const-min_,mult,g0)
            ys0 = ansatz(xs, *p0)

            # to_plot = r in plot
            to_plot = r % 20 == 0 or ii == 0

            try:
                p1, pcov = scipy.optimize.curve_fit(ansatz, xs, ys, p0)
                    # ,bounds=([0,0,0],[np.inf,np.inf,1000]))
                perr = np.sqrt(np.diag(pcov))
                ys1 = ansatz(xs, *p1)

                rs.append(r)
                consts.append(p1[0])
                mult1s.append(p1[1])
                if len(p1) > 2:
                    weight1s.append(p1[2])
                if len(p1) == 4:
                    as_.append(p1[3])
                if len(p1) == 5:
                    mult2s.append(p1[3])
                    weight2s.append(p1[4])
                errs.append(perr/p1)

                # print(r,0,0,list(ys1))
                break

            except Exception as e:
                ys1 = p1 = perr = None
                if r > lim_:
                    print('error', r, e)

        if to_plot and show_fit_plots:
            # print('%d: ~%r ... %r±%r' % (r, p0, p1, perr))
            fig, ax = plt.subplots()
            ax.plot(xs, ys)
            ax.plot(xs, ys0, ':')
            if ys1 is not None:
                ax.plot(xs, ys1)
            ax.grid()
            ax.set_xlim(-r,r)
            ax.set_title(repr(p0)+"\n"+repr(p1))
            plt.show(fig)
            plt.close(fig)

    fig, ax = plt.subplots()
    ax.semilogy(rs, errs)
    ax.set_title('errors')
    fig, ax = plt.subplots()
    ax.semilogy(rs, consts)
    ax.set_title('consts')
    fig, ax = plt.subplots()
    ax.semilogy(rs, mult1s)
    ax.set_title('mult1s')
    # fig, ax = plt.subplots()
    # ax.semilogy(rs[:lim], [m-c for c,m in zip(consts,mult1s)][:lim])
    # ax.set_title('consts\'')
    if weight1s:
        fig, ax = plt.subplots()
        ax.plot(rs, weight1s)
        ax.set_title('weight1s')
        # fig, ax = plt.subplots()
        # ax.plot(rs, [(np.exp(w)-1 if w > 0 else w) for w in weight1s])
    if as_:
        fig, ax = plt.subplots()
        ax.plot(rs, as_)
        ax.set_title('as_')
    if mult2s:
        fig, ax = plt.subplots()
        ax.semilogy(rs, mult2s)
        ax.set_title('mult2s')
        fig, ax = plt.subplots()
        ax.plot(rs, weight2s)
        ax.set_title('weight2s')
    plt.show()
    plt.close('all')
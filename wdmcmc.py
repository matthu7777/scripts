#!/usr/bin/python

# builtins
import sys
import math, os, random
import os.path

# 3rd party
import numpy as np
import scipy.linalg
from ppgplot import *
import emcee

# mine
import trm.subs as subs
import trm.subs.input as inp
import trm.mcmc as mcmc
from   trm.subs import Fname

usage = """usage: %s data model init prior chain [append [lmarq nscale scale cburn cmin nstore cstart nstart]] nrate ntrial

Program to run MCMC on lcurve models. This runs by random jumps in model space
which are chosen or not depending upon a posterior probability evaluated
through a combination of chi**2 (the prob(data given model)) and a prior model
that is defined by user-supplied code.

This is essentially a management task with all the hard work hived off to
lroche. The great difficulty is to jump sensibly and maintain a decent success
rate. This script tries to handle this, but only in a fairly generic way and
you should always be aware of the possible difficulties that can occur. See at
the end of this for a discussion of these.


Arguments:

data    -- file of lightcurve data for input to levmarq, lroche, etc. All the
           various files have standard 3 letter suffixes which will be added
           on, e.g. '.dat' for data so that you can specify the same name
           (without the suffix) for most of them without fear of conflict.

model   -- light curve model which fits data. Extension '.mod'

init    -- file of ranges of each variable used to initialise the
           chain. Variables are picked uniformly around starting values in the
           model file

           Example format: 

           q = 0.01 
           r1 = 0.01 

           where the numbers are the half ranges to perturb each parameter
           by. Any variable not specified will be set to its value in
           'model'. Careful use of this file can locate hard-to-find minima in
           chi**2 space but also implies an initial 'burn-in' period for the
           chain that can be long or even infinite if it gets stuck in some
           high hidden valley somewhere. I have found this to be more than
           likely in models I was running, but it will depend greatly on which
           parameters you are varying.

prior --   python file which implements the user-defined prior. This file, if
           specified, must contain a function called 'prior' which returns
           ln(prior probability). This makes a gaussian prior easy to
           implement as -((v-m)/sigma)**2/2.  The function is fed the parameter
           names and values via a dictionary, its only argument. Here is an
           example enforcing an inclination close to 85 degrees:

           def prior(p):
              return -((p['iangle']-85)/0.1)**2/2.

           A more complex one could enforce a near M-S mass-radius relation
           for example. Note that normalisation constants (sqrt(2\pi) and the
           like) are irrelevant. Here is a prior enforcing preventing
           epow_spot from getting too large and simultaneously keeping the
           disk exponent in check

           def prior(p):
              lp = 0.
              if epow_spot > 3.:
                 lp += (p['epow_spot']-3.)**2
              lp += ((p['text_disc']-1.)/0.5)**2
              return -lp/2.

           An invalid file name will result in no prior being applied.

chain --   output file for storage of results. Extension '.log'. If the file
           already exists, then the next parameter gives you the option to use
           it to define the chain and to append new results to it.

append --  If the file specified by 'chain' already exists, it can be used as
           the basis for continuing a run. In this case the information
           necessary to define the iterations is read from the file which is
           opened for appending. If the file does not exist, this parameter
           will not be prompted for and will be set = False

If append = True:

lmarq --   levmarq-like output specifying the covariances for the jump
           distribution. Extension '.lmq'. This can either come from levmarq
           itself or the script genlmq.py. See at the end for a detailed
           discussion of this. The script uses the covariances to generate
           jumps using an equivalent multi-variate normal distribution.

nscale --  the jumps are scaled by a factor which is updated after every nscale
           successes aiming at a success rate of 0.25. The scale factor is
           multiplied by 4*nscale/ntrial where ntrial is the number of trials
           taken to obtain the nscale successes.

scale --   the jumps are scaled by this factor to allow you to vary the
           acceptance rate. This is the initial value.

cmin --    minimum chi**2. This number is used to re-scale the uncertainties to
           get chi**2 = 1 per dof before adding in the prior probability
           factor 'lp'. If you really believe your uncertainties, set this
           equal to the weighted number of degrees of freedom reported by
           'lroche' (i.e. 'the value of 'wnok'). Quite often cmin is larger
           than wnok because of correlated noise (e.g. flickering) in the
           light curve, and one cannot expect the crude re-scaling to be the
           only fix. Just do your best under such circumstances.

nstore  -- how often to store results to chain

In all cases:

nrate   -- how often to print out the acceptance rate. Default = 1000

ntrial  -- maximum number of trials

Jump distributions
==================

The most difficult and important part here is the development of a good jump
distribution. There is no single solution here.  First you might want to
generate an artificial (diagonal) jump distribution with genlmq.py if you
don't have one already. Do your best to supply reasonable estimates for the
sigmas here because if you are grossly off in their relative magnitudes, the
calculations will be slowed down. Then carry out a shortish run to generate a
chain from which to compute a better approximation.  Apply covar.py to this to
generate a better jump distribution and then re-run. Once you are happy with
the jump distribution, go for a 'production' run where the jump distribution
is held fixed. The big difficulty is knowing when to stop developing the jump
distribution and some of the pitfalls this entails. The major difficulties are
caused by strongly correlated parameters. These are difficult to start with
when one uses uncorrelated jumps which have to be very small not to jump
uphill. Slow random-walk, red-noise type chains with long correlation lengths
(as opposed to the ideal white noise appearance) are the consequence of
this. Note that it only may appear like this in some, not all parameters. Use
pall.py, series.py, twod.py, corr.py to diagnose these. Once you have got an
lmq file reflecting such correlations another problem can emerge: you can end
up with overly strong correlations as the result of exploring too limited a
part of a narrow winding path in parameter space. This can make it very hard
ever to reach regions of parameter space not on your original track. In
principle the gaussian jump distribution can reach anywhere, and detailed
balance ensures that the chain will represent the correct distribution, but it
is entirely possible to set it up so that in practice this will almost never
happen. The other side of the coin is that if a chain does reach such a
region, it can get stuck there for extended periods. The result is a chain
that has long sections of what looks close to white noise and is well mixed,
interspersed with jumps to fairly different regions. The trouble is, such
jumps may be so rare that you don't even realise that they can occur giving a
spuriously rosy view of the parameter uncertainties. Two ways to combat this
are (1) to start multiple chains off at different places, and (2) to temper
the correlations using the script modlmq.py which pushes the distribution
towards an uncorrelated one. Be very wary: it is only too easy to come a
cropper here. It is very easy to lose patience and assume that your chain has
converged. Try, if you can, to speed your models so that you can run long
chains to test this. 100,000 or more light curve computations may be
needed. Blunting the correlations can have a price of course, as it moves you
back towards the uncorrelated distribution where uphill jumps are common ad
motion in parameter space diffusive in character, but it does raise the
chances of an "unexpected" step.

"""

# print help
if len(sys.argv) == 2 and sys.argv[1] == '-h':
    print usage % (sys.argv[0][sys.argv[0].rfind('/')+1:],)
    exit(0)

# initiate global variable for white dwarf contribution (MG)
# necessary because wdwarf cannot be passed out of Lpprob()
#global_wd = 0

# generate arguments
inpt = inp.Input('PYTHON_MCMC_ENV', '.pymcmc', sys.argv)

# register parameters
inpt.register('data',    inp.Input.LOCAL, inp.Input.PROMPT)
inpt.register('model',   inp.Input.LOCAL, inp.Input.PROMPT)
inpt.register('init',    inp.Input.LOCAL, inp.Input.PROMPT)
inpt.register('chain',   inp.Input.GLOBAL,inp.Input.PROMPT)
inpt.register('append',  inp.Input.GLOBAL,inp.Input.PROMPT)
inpt.register('prior',   inp.Input.LOCAL, inp.Input.PROMPT)
inpt.register('method',  inp.Input.LOCAL, inp.Input.PROMPT)
inpt.register('nwalker', inp.Input.LOCAL, inp.Input.PROMPT)
inpt.register('stretch', inp.Input.LOCAL, inp.Input.PROMPT)
inpt.register('lmarq',   inp.Input.LOCAL, inp.Input.PROMPT)
inpt.register('nscale',  inp.Input.LOCAL, inp.Input.PROMPT)
inpt.register('scale',   inp.Input.LOCAL, inp.Input.PROMPT)
inpt.register('cmin',    inp.Input.LOCAL, inp.Input.PROMPT)
inpt.register('nstore',  inp.Input.LOCAL, inp.Input.PROMPT)
inpt.register('nrate',   inp.Input.LOCAL, inp.Input.HIDE)
inpt.register('ntrial',  inp.Input.LOCAL, inp.Input.PROMPT)

# Get inputs
data      = inpt.get_value('data',   'file with light curve data', Fname('ippeg', '.dat'))
model     = inpt.get_value('model',  'file with converged light curve model for this data', Fname('ippeg', '.mod'))
init      = inpt.get_value('init',   'file with initial ranges for variable parameters', Fname('ippeg', '.ini'))
fptr  = open(init)
idict = {}
for line in fptr:
    if not line.startswith('#'):
        eq = line.find('=')
        if eq > -1:
            par   = line[:eq].strip()
            idict[par] = float(line[eq+1:].strip())
fptr.close()

# Prior constraints can defined using a piece of user-defined code
try:
    prior = inpt.get_value('prior',  'python code defining prior constraints on parameters', Fname('ippeg', '.py'))  
    sys.path.insert(0,'.')
    Prior = __import__(prior[:-3])
    prior = Prior.prior

except inp.InputError:
    prior = None
    print 'No prior will be applied.'

# create the model
lcurve = mcmc.Lcurve(model, prior)

# get the name of the log file.
chain    = inpt.get_value('chain',  'file for storage of mcmc output', Fname('ippeg', '.log', Fname.NEW))
append   = chain.exists() and inpt.get_value('append',  'append output to this file (otherwise overwrite)?', True)

if not append:
    method  = inpt.get_value('method',   'method to use, "m" for Metroplis or "a" for affine', \
                             'm', lvals=['a','m'])

    if method == 'a':
        nwalker = inpt.get_value('nwalker', 'number of walkers',100, max(2,2*lcurve.nvar()), multipleof=2)
        stretch = inpt.get_value('stretch', 'stretch factor', 2., 1.)
    else:
        nwalker = 1
        stretch = 0.
        nsub    = 1

    # need to get several other parameters if we are starting from scratch
    lmarq   = inpt.get_value('lmarq',  'file with levmarq output for covariance data', Fname('ippeg', '.lmq', Fname.OLD))

    jump    = mcmc.Jump(lmarq)
    # cut down covariances to match model variables
    jump.trim(lcurve.var().keys())
    scale   = inpt.get_value('scale', 'initial scale factor', 0.1, 1.e-10)
    jump   *= scale

    cmin    = inpt.get_value('cmin',   'minimum chi**2', 100., 0.)
    nstore  = inpt.get_value('nstore', 'how often to store results', 10, 1)

inpt.set_default('nrate', 1000)
nrate   = inpt.get_value('nrate', 'how often to print acceptance rate', 1000, 1)

# always need the following parameters
ntrial  = inpt.get_value('ntrial', 'number of trials', 1000, 1)

nrep = 50
print 'Will report progress 1 in every ' + str(nrep) + ' trials and each time an improved model is produced'

if append:

    # in this case we use the last line to set the model for which we have to 
    # read the names from the starting lines
    clog = mcmc.Fchain(chain)
    print 'Read chain from',chain
    lvars = clog.vals[-1,:-3]
    if len(lvars) == lcurve.nvar():
        for (par,val) in zip(clog.names[:-3],lvars):
            lcurve[par][0] = val
    else:
        print 'Number of values of last line of ' + chain + ' = ' + str(len(lvars)) + \
            ' which differs from the number of variables read from the model file = ' + str(lcurve.nvar())
        exit(1)

    cold,wnok,wdold = lcurve.chisqlp(data)
    lpold = lcurve.lnprior()
    print 'Initial chi**2,lp,wnok,wdold =',cold,lpold,wnok,wdold
    jump    = clog.jump
    cmin    = clog.chmin
    nstore  = clog.nstore
    method  = clog.method
    nwalker = clog.nwalker
    stretch = clog.stretch
    if len(clog) % nwalker != 0:
        print 'The number of lines =',len(clog.vals),'was not a multiple of nwalker =',nwalker
        exit(1)

else:

    # starting without a log file; need first to
    # find an acceptable starting model
    ntry = 0
    ok = False
    cvar = lcurve.var()
    while not ok and ntry < 1000:
        for (par,hrange) in idict.iteritems():
            if par in cvar:
                lcurve[par][0] = random.uniform(cvar[par]-hrange,cvar[par]+hrange)
        ok = lcurve.ok()
        if ok:
            cold,wnok,wdold = lcurve.chisqlp(data)
            ok = cold is not None
        ntry += 1

    if not ok:
        print 'Tried',ntry,'initial models, none of which worked; giving up.'
        exit(1)

    lpold = lcurve.lnprior()
    print 'OK initial model model found after',ntry,'trials.'
    print 'Initial chi**2,lp,wnok,wdold =',cold,lpold,wnok,wdold
    model = dict([(key, lcurve[key]) for key in lcurve])
#    clog  = mcmc.Fchain(model, nstore, method, jump, nwalker, stretch, ['chisq', 'lp', 'wdwarf', 'pprob'], cmin, chain)
    clog  = mcmc.Fchain(model, nstore, method, jump, nwalker, stretch, ['chisq', 'wdwarf', 'lp', 'pprob'], cmin, chain)

vvar = clog.vpars()

def lnpost(p, lcurve, data, cmin):
    """
    Computes ln(posterior prob) given a trial model, an Lcurve model,
    a data file name and a minimum chi**2 for error bar re-scaling.
    Returns: ln(post), ln(prior), chisq. If anything fails, it comes
    back with -Inf,0,-Inf
    """

    lcurve.set(p)
    if lcurve.ok():
        lnprior = lcurve.lnprior()
        chisq, wnok, wdwarf = lcurve.chisqlp(data)

        # save wdwarf to a global variable so that it can be passed out of Lpprob()
        # (MG)
        #global global_wd
        #global_wd = wdwarf
        #print wdwarf, global_wd ##

        if chisq is not None:
            # scale chisq according to the preset cmin value
            lnpost = -wnok/cmin*chisq/2.+lnprior
            return lnpost, lnprior, chisq
        else:
            return (float('-Inf'),0.,float('Inf'))
    else:
        return (float('-Inf'),0.,float('Inf'))

class Lpprob(object):
    """This returns the ln of posterior probability as required by emcee in
    'function object' form given vector x of variable values. It is assumed
    that the variables in x run in the same order as in vvar, the names of the
    variables. 'lcurve' represents the model and prior 'data' is the name of
    the data file. 'cmin' is the minimum chi-squared predicted which is used
    to re-scale error bars. Any failures are converted into very low
    probabilities.

    If method == 'a' then a straight set of values is expected in the call to
    this object, else a dictionary

    """

    def __init__(self, lcurve, vvar, data, cmin, method):
        self.lcurve = lcurve
        self.vvar   = vvar
        self.data   = data
        self.cmin   = cmin
        self.method = method

    def __call__(self, x):
        if method == 'a':
            p  = dict(zip(self.vvar, x))
        else:
            p = x

        lnpo, lnpri, chisq = lnpost(p, self.lcurve, self.data, self.cmin)
        return lnpo

# OK get going
mlpost   = cold*wnok/cmin + lpold
nt       = 0
nsuccess = 0
nsuc     = 0
ntry     = 0

# create the posterior probability function
lnprob = Lpprob(lcurve, vvar, data, cmin, method)

if method == 'a':

    if clog.vals is None or not len(clog.vals):
        # Generate the walkers. We use the jump distribution
        # to make a cloud of (ok) vectors
        walkers = []
        pinit   = lcurve.var()
        for i in xrange(nwalker):
            ok = False
            while not ok:
                ptest = jump.trial(pinit)
                lcurve.set(ptest)
                ok = lcurve.ok()
            walkers.append(np.array([ptest[vnam] for vnam in vvar]))
    else:
        # Generate walkers from the file just read
        walkers = clog.vals[-nwalker:,:-3]

    sampler = emcee.EnsembleSampler(nwalker, len(vvar), lnprob, a=stretch)

    lnpmax = -1.e30

    # Do the hard work
    lnps = None
    rs   = None
    ntot = 0
    for i in xrange(ntrial/nstore):

        # go through in steps to reduce memory requirement
        mrate = 0.
        walkers, lnps, rs = sampler.run_mcmc(walkers, nstore, rs, lnps)

        ntot  += nwalker * nstore
        print 'Group',i+1,'acceptance rate =',sampler.acceptance_fraction.mean()
        sys.stdout.flush()

        # check to see if we have improved the model
        if lnps.max() > lnpmax:
            lnpmax = lnps.max()
            print 'New best model has ln(post) =',lnpmax

        sampler.reset()

        # Write the current position to a file, one line per walker, computing
        # the chisq and lnpriors
        for walker, lnp  in zip(walkers, lnps):
            p     = dict(zip(vvar, walker))
            lcurve.set(p)
            lnpri = lcurve.lnprior()
            chisq = 2.*cmin*(lnpri - lnp)/wnok

            # read white dwarf contribution from global variable (MG)
            #wdwarf = global_wd
            temp1, temp2, wdwarf = lcurve.chisqlp(data)

            clog.add_line(np.concatenate((walker, [chisq, wdwarf, lnpri, lnp])))

else:

    # Main long loop
    while nt < ntrial:
        nt   += 1
        ntry += 1

        # main computation part
        success,cold,lpold,wdold = lcurve.trial(
            jump,data,cold,lpold,cmin,None,wdold)
        if success:
            nsuc += 1
            nsuccess += 1

        if nt % nrep == 0:
            print 'Chi**2,lp,mean success rate after trial',nt,\
                '=',cold,lpold,nsuccess/float(nt)

        if ntry == nrate:
            clog.add_info('# success rate after trial ' + str(nt)
                          + ' = ' + str(nsuc/float(ntry)) + '\n')
            nsuc = 0
            ntry = 0

        mlpos = cold*wnok/cmin + lpold

        # save model to file
        if nt % nstore == 0:
            clog.add_line(lcurve.var().values() + [cold,lpold,wdold,mlpos])


#!/usr/bin/python

from pyevolve import G1DList
from pyevolve import GSimpleGA
from pyevolve import Selectors
from pyevolve import Mutators
from pyevolve import Crossovers
from pyevolve import Initializators
from pyevolve import GAllele
from pyevolve import Consts
from pyevolve import DBAdapters
from pyevolve import Scaling
from numpy import genfromtxt, array, arange, linalg, exp, diff, concatenate, round_
from numpy import where, round as rnd, savetxt, loadtxt, std, mean, linspace, ceil, min, max
from scipy.stats import norm
from commands import getstatusoutput as gso
from string import ascii_letters, digits
from random import sample
from glob import glob
import pylab as plt
import os
import sys
from initial_mpl import init_plotting_isi
from PyHyp71 import run_hyp71

"""

Script for calculating 1D-velocity model using GA algorithm.

Note:

Using Hypo71 as an objective function still is in progress.

ChangeLogs:

07-Aug-2017 > Initial.
23-Aug-2017 > Set depth to 0/.5 decimal values.
02-Oct-217  > Fixed some issues in eval_func(), removed "23-Aug-2017" change, fixed issues for plotting.

"""
    
#________________WRITE MODELS INTO DB

def model_writer(file_name, model):

    with open(file_name,'a') as f:

        savetxt(f, model, fmt='%5.2f', newline=" ")
        f.write('\n')
        
#________________SET INITIAL PARAMETERS

inp_dic = {}

try:

    inp = genfromtxt('par.dat', delimiter='=', dtype=None)

except IOError:

    with open('par.dat', 'w') as f:

        f.write("""#
#
# INPUT FILE PARAMETER FOR 'PyGAVel' v0.1.
#
###################
#
SYNTHETIC_F        = True # Synthetic test enalbe/disable.
REAL_MODEL_V       = 4.50, 5.20, 5.50, 5.80, 6.50
REAL_MODEL_D       = 0.00, 4.00, 8.00,12.00,18.00
REAL_MODEL_R       = 1.83, 1.80, 1.75, 1.70, 1.70
#
MODEL_VEL_MIN      = 4.00, 4.00, 4.00, 4.00, 4.00, 4.00, 4.00, 4.00, 4.00, 4.00
MODEL_VEL_MAX      = 7.00, 7.00, 7.00, 7.00, 7.00, 7.00, 7.00, 7.00, 7.00, 7.00
MODEL_DEP_MIN      = 0.00, 2.00, 4.00, 6.00, 8.00,10.00,12.00,14.00,16.00,18.00
MODEL_DEP_MAX      = 0.00, 2.00, 4.00, 6.00, 8.00,10.00,12.00,14.00,16.00,18.00
MODEL_VPS_MIN      = 1.83, 1.83, 1.80, 1.80, 1.75, 1.75, 1.70, 1.70, 1.70, 1.70
MODEL_VPS_MAX      = 1.83, 1.83, 1.80, 1.80, 1.75, 1.75, 1.70, 1.70, 1.70, 1.70
#
GENSIZE            = 2
POPSIZE            = 2
MUT_R              = 0.10
MUT_R_V            = 0.00 # Mutation decreasing rate. Different "MUT_R" for each generation. MUT_R=exp(-MUT_R_V*range(GENSIZE)). "0" = Disable.
CRO_R              = 0.90
NUM_ELTSM          = 1
SEED               = 101
NUM_CPU            = 4    # Set 0 to disable it.
OBJ_FUNC           = 2    # 1: Hypo71, 2: Hypoelipse. Note to use POS_GRAD_V=True if Hypo71 is selected.
POS_GRAD_V         = True
DBASE              = noisy_vpvs_1
RESET_DB           = True
GA_SELECTOR        = T    # T:Turnoment, R:RouletWhell, N:Rank, M: T(75%)+ R(20%) and N(5%).
#
PLT_VP_RNG         = 3.00, 8.00
PLT_DP_RNG         = 0.00, 25.0
PLT_R_RNG          = 1.65, 1.85
""")
        print '\n+++ No "par.dat" file was found.'
        print '+++ The parameter file "par.dat" has been created. Modify it and run again.'
        sys.exit(0)

for line in inp:

    try:

        inp_dic[line[0].strip()] = float(line[1])

    except ValueError:

        if len(line[1].split(',')) > 1:

            inp_dic[line[0].strip()] = array(line[1].split(','), dtype=float).tolist()

        else:

            inp_dic[line[0].strip()] = line[1].strip()

if eval(inp_dic['SYNTHETIC_F']):
 
     velocity_real = inp_dic['REAL_MODEL_V']
     depth_real    = inp_dic['REAL_MODEL_D']
     vpvs_real     = inp_dic['REAL_MODEL_R']

velocity_min   = inp_dic['MODEL_VEL_MIN']
velocity_max   = inp_dic['MODEL_VEL_MAX']
depth_min      = inp_dic['MODEL_DEP_MIN']
depth_max      = inp_dic['MODEL_DEP_MAX']
vpvs_min       = inp_dic['MODEL_VPS_MIN']
vpvs_max       = inp_dic['MODEL_VPS_MAX']
generationSize = inp_dic['GENSIZE']
populationSize = inp_dic['POPSIZE']
pMutation      = inp_dic['MUT_R']
mut_rv         = inp_dic['MUT_R_V']
pCrossover     = inp_dic['CRO_R']
num_etlsm      = inp_dic['NUM_ELTSM']
pos_grad_v     = inp_dic['POS_GRAD_V']
dbase_name     = inp_dic['DBASE']
resetDB        = inp_dic['RESET_DB']
num_cpu        = inp_dic['NUM_CPU']
obj_func       = inp_dic['OBJ_FUNC']
seed           = inp_dic['SEED']
ga_selector    = inp_dic['GA_SELECTOR']
random_pool    = ascii_letters+digits
mut_list       = exp(-arange(generationSize)*mut_rv)
models_file    = open('models.dat','w')

mut_list[mut_list<pMutation] = pMutation

#________________OBJECTIVE FUNCTION

def eval_func(chromosome):

    ID = ''.join(sample(random_pool,6))

    #__________SORT VELOCITIES IF POSITIVE GRADIENT IS REQUESTED

    if eval(pos_grad_v):

        chromosome[:len(velocity_min)] = sorted(chromosome[:len(velocity_min)])
        chromosome[:len(velocity_min)] = concatenate(([0],0<abs(where(diff(chromosome[:len(velocity_min)])<.05,1,0))))*.05+chromosome[:len(velocity_min)]

    #__________WRITE INDIVIDUALS INTO DB

    model_writer(models_file.name, chromosome.genomeList)

    #__________MAKE AN ABRUBT CHANGE IN MUTATION RATE IF REQUIRED

    if mut_rv > 0.0: ga.pMutation = mut_list[ga.getCurrentGeneration()]

    #__________START HYPO71

    if obj_func == 1:

        hyp71_pha = 'norhyp.out'
        nord_sta  = 'STATION0.HYP'
        vel_model = array([chromosome[:len(velocity_min)],chromosome[len(velocity_min):2*len(velocity_min)]]).T
        score     = run_hyp71(hyp71_pha_file=hyp71_pha, nordic_sta_file=nord_sta, vel_model=vel_model, run_id=ID)
        
    #__________START HYPOELLIPSE

    if obj_func == 2:
     
        with open(ID+'.prm','w') as f:

            for v,d,r in zip(chromosome[:len(velocity_min)],
                             chromosome[len(velocity_min):2*len(velocity_min)],
                             chromosome[2*len(velocity_min):]):

                f.write('VELOCITY            %5.2f %5.2f %4.2f\n'%(v,d,r))

        cmd = 'fwd_problem.sh '+ID+' > /dev/null'
        os.system(cmd)
        score = genfromtxt(ID+'misfit.val')

        for _ in glob(os.path.join(ID+'*')): os.remove(_)

    #__________RETURN SCORE

    return score

#________________RUN GA

def run_ga():
    
    global ga

    #___________________Genome instance
    # 

    setOfAlleles = GAllele.GAlleleList()
  
    pars_min = velocity_min + depth_min + vpvs_min
    pars_max = velocity_max + depth_max + vpvs_max
    num_pars = len(pars_min)

    for (vmin, vmax) in zip(pars_min, pars_max):

        tmp = GAllele.GAlleleRange(vmin, vmax, real=True)
        setOfAlleles.add(tmp)

    genome = G1DList.G1DList(num_pars)
    genome.setParams(allele=setOfAlleles)

    genome.initializator.set(Initializators.G1DListInitializatorAllele)
    genome.mutator.set(Mutators.G1DListMutatorAllele)
    genome.crossover.set(Crossovers.G1DListCrossoverUniform)

    #___________________The evaluator function (objective function)
    #

    genome.evaluator.set(eval_func)

    #___________________Genetic Algorithm Instance
    #

    ga = GSimpleGA.GSimpleGA(genome, seed=int(seed))

    if num_cpu: ga.setMultiProcessing(True, True, int(num_cpu))

    if ga_selector == 'T': ga.selector.set(Selectors.GTournamentSelector)
    if ga_selector == 'R': ga.selector.set(Selectors.GRouletteWheel)
    if ga_selector == 'N': ga.selector.set(Selectors.GRankSelector)

    if ga_selector == 'M':

        ga.selector.setRandomApply(True)
        ga.selector.set(Selectors.GTournamentSelector,0.75)
        ga.selector.add(Selectors.GRouletteWheel,0.20)
        ga.selector.add(Selectors.GRankSelector)

    ga.setMinimax(Consts.minimaxType["minimize"])
    ga.setGenerations(int(generationSize))
    ga.setPopulationSize(int(populationSize))
    ga.setCrossoverRate(pCrossover)
    ga.setMutationRate(pMutation)
    ga.setElitism(True) 
    ga.setElitismReplacement(int(num_etlsm))

    #___________________Sets the DB Adapter
    #

    sqlite_adapter = DBAdapters.DBSQLite(identify=dbase_name, resetDB=eval(resetDB))
    ga.setDBAdapter(sqlite_adapter)

    #___________________Do the evolution
    #

    ga.evolve(freq_stats=5)

    #___________________Print Best individual
    #

    best    = ga.bestIndividual()
    best_rs = best.getRawScore()
    best_v  = best.genomeList[:len(velocity_min)]
    best_d  = best.genomeList[len(velocity_min):2*len(velocity_min)]
    best_r  = best.genomeList[2*len(velocity_min):]

    print ''
    print '+++ Best Raw Score =',best_rs
    print '+++ FinalModel :'
    print '   +++ Velocity :',rnd(best_v,2)
    print '   +++ Depth    :',rnd(best_d,2)
    print '   +++ VpVs     :',rnd(best_r,2)

    return best, best_rs, best_v, best_d, best_r

#________________WRITE FINAL RESULT

def write_res(dbase_name, flag, best_rs, best_v, best_d, best_r):

    if flag: res = open('result.dat','w')
    else: res = open('result.dat','a')

    res.write('project_name:%s; best raw score=%7.4f\n'%(dbase_name, best_rs))

    for v,d,r in zip(best_v, best_d, best_r):

        res.write('%5.2f %5.2f %5.2f\n'%(v,d,r))

    res.close()

#________________PLOT RESULTS

def plot(best, best_rs, best_v, best_d, best_r):

    init_plotting_isi(16,8)

    #___________________Plot final results
    #

    if eval(inp_dic['SYNTHETIC_F']):
        
        vel_list  = [velocity_min, velocity_max, best_v, velocity_real]
        dep_list  = [depth_min, depth_max, best_d, depth_real]
        vpvs_list = [vpvs_min, vpvs_max, best_r, vpvs_real]
        colors    = ['k','k','b','r']
        labels    = ['Min', 'Max', 'Best', 'Real']

    else:

        vel_list  = [velocity_min, velocity_max, best_v,]
        dep_list  = [depth_min, depth_max, best_d]
        vpvs_list = [vpvs_min, vpvs_max, best_r]
        colors    = ['r','k','g']
        labels    = ['Min', 'Max', 'Best']

    ax = plt.subplot(121)
    [i.set_linewidth(0.6) for i in ax.spines.itervalues()]

    for v,d,c,l in zip(vel_list, dep_list, colors, labels):
     
        xs = []
        ys = []

        x = array(v)
        y = array(d)

        for i,j in zip(x,y):

            xs.append(i)
            xs.append(i)
            ys.append(j)
            ys.append(j)

        xs.pop(-1)
        ys.pop(0)
        xs.append(xs[-1])
        ys.append(max(inp_dic['PLT_DP_RNG']))

        if l == 'Min': ax.plot(array(xs),-array(ys), linewidth=1.5, color='k', linestyle='--', label=l)
        elif l == 'Max': ax.plot(array(xs),-array(ys), linewidth=1.5, color='k', linestyle='-.', label=l)
        else: ax.plot(array(xs),-array(ys), linewidth=1.5, color=c, linestyle='-', label=l)
   
    ax.set_xlabel('Velocity [km/s]')
    ax.set_ylabel('Depth [km]')
    ax.set_xlim(inp_dic['PLT_VP_RNG'])
    ax.set_ylim(-array(inp_dic['PLT_DP_RNG'])[::-1])
    ax.locator_params(axis = 'x', nbins = 6)
    ax.locator_params(axis = 'y', nbins = 6)
    ax.grid()
    ax.legend(loc=3)

    ax = plt.subplot(122)
    [i.set_linewidth(0.6) for i in ax.spines.itervalues()]

    for r,d,c,l in zip(vpvs_list, dep_list, colors, labels):
     
        xs = []
        ys = []

        x = array(r)
        y = array(d)

        for i,j in zip(x,y):
             
            xs.append(i)
            xs.append(i)
            ys.append(j)
            ys.append(j)

        xs.pop(-1)
        ys.pop(0)
        xs.append(xs[-1])
        ys.append(max(inp_dic['PLT_DP_RNG']))

        if l == 'Min': ax.plot(array(xs),-array(ys), linewidth=1.5, color='k', linestyle='--', label=l)
        elif l == 'Max': ax.plot(array(xs),-array(ys), linewidth=1.5, color='k', linestyle='-.', label=l)
        else: ax.plot(array(xs),-array(ys), linewidth=1.5, color=c, linestyle='-', label=l)

    ax.set_xlabel('VpVs [km/s]')
    ax.set_ylabel('Depth [km]')
    ax.set_xlim(inp_dic['PLT_R_RNG'])
    ax.set_ylim(-array(inp_dic['PLT_DP_RNG'])[::-1])
    ax.locator_params(axis = 'x', nbins = 6)
    ax.locator_params(axis = 'y', nbins = 6)
    ax.grid()
    ax.legend(loc=2)

    plt.tight_layout()
    plt.savefig(dbase_name+'.tiff',dpi=300)
    plt.close()

    #__________StdDev (V,D,R)
    
    models    = loadtxt('models.dat')
    model_std = std(models,axis=0)[0]
    model_men = mean(models,axis=0)[0]
    best      = array([best_v,best_d,best_r]).flatten()
    tot_ax    = ceil(models.shape[1]/4.)

    init_plotting_isi(16,16)
    plt.rcParams['xtick.labelsize'] = 6
    plt.rcParams['ytick.labelsize'] = 6
    plt.rcParams['axes.labelsize']  = 6
    
    for par in range(models.shape[1]):

        ax = plt.subplot(tot_ax,4,par+1)
        [i.set_linewidth(0.6) for i in ax.spines.itervalues()]

        ax.text(0.50, 1.18, 'Parameter-%d'%(par+1), fontsize=6,
                transform=ax.transAxes,ha='center', va='top')

        mu    = mean(models[:,par])
        sig   = std(models[:,par])
        x     = linspace(mu-3*sig,mu+3*sig, 100)
        data  = models[:,par]
        a,b,c = plt.hist(data, 50, color='r', linewidth=0, alpha=.6)
        plt.vlines(best[par],0,max(a),color='g',zorder=20,linewidth=2)
        plt.vlines(mu,0,max(a),color='b',zorder=20,linewidth=2)
        plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
        plt.locator_params(axis = 'x', nbins = 6)
        plt.locator_params(axis = 'y', nbins = 4)

    plt.tight_layout(True)
    plt.savefig('models_stat.tiff',dpi=300)
    plt.close()

#________________DO

import warnings
warnings.filterwarnings("ignore")

if mut_rv > 0.0:

    fig = plt.figure()
    fig.set_tight_layout(True)

    ax1 = plt.subplot(111)
    ax1.plot(mut_list, 'ro-')
    ax1.set_xlabel('Generations #')
    ax1.set_ylabel('Mutation Rate')
    ax1.locator_params(axis = 'x', nbins = 6)
    ax1.locator_params(axis = 'y', nbins = 6)
    ax1.grid()
    plt.savefig(dbase_name+'_mut_r.tiff',dpi=300)
    plt.close()

best, best_rs, best_v, best_d, best_r = run_ga()
write_res(dbase_name, eval(resetDB), best_rs, best_v, best_d, best_r)
plot(best, best_rs, best_v, best_d, best_r)

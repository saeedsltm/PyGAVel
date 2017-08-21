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
from numpy import genfromtxt, array, arange, linalg, exp, diff, concatenate
from numpy import where, round as rnd, savetxt, loadtxt, std, mean, linspace, ceil
from scipy.stats import norm
from commands import getstatusoutput as gso
from string import ascii_letters, digits
from random import sample
from glob import glob
import pylab as plt
import os
import sys
from initial_mpl import init_plotting
from PyHypo71 import run_hyp71

"""
Script for calculating 1D-velocity model using GA algorythm.

Note:

Using Hypo71 as an objective function still is in progress.

ChangeLogs:

07-Aug-2017 > Initial.

"""

#________________PLOT DATA

def hist_depth():

    if not os.path.exists('hypoel.sum'):

        print '\n\n+++ No "hypoel.sum" was found!'
        sys.exit(0)

    dep = []

    with open('hypoel.sum') as f:

        for l in f:

            dep.append(-float(l[31:36])*1e-2)

    init_plotting()

    ax   = plt.subplot(111)
    bins = arange(min(dep), max(dep) + 1, 1)
    ax.hist(dep,bins=bins,orientation="horizontal",color='r',alpha=.7)
    ax.set_xlabel('Number of Events [#]')
    ax.set_ylabel('Depth [km]')
    ax.grid(True)
    plt.show()
    plt.close()

##hist_depth()
##
##ans = raw_input('\n+++ Press "Enter" to continue, "q" for quit.\n')
##
##if ans == 'q':
##
##    sys.exit(0)
    
#________________WRITE MODELS INTO DB

def model_writer(file_name, model):

    with open(file_name,'a') as f:

        savetxt(f, model, fmt='%5.2f', newline=" ")
        f.write('\n')
        
#________________SET INITIAL PARAMETERS

inp_dic = {}
inp     = genfromtxt('par.dat', delimiter='=', dtype=None)

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

    #__________SORT VELOCITIES IF REQUIRED

    if eval(pos_grad_v):

        for _ in xrange(ga.getPopulation().popSize):

            ga.getPopulation()[_].genomeList[:len(velocity_min)] = sorted(ga.getPopulation()[_].genomeList[:len(velocity_min)])

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

        cmd = './fwd_problem.sh '+ID+' > /dev/null'
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
    print '+++ Best Raw Score=',best_rs
    print '+++ FinalModel:\n'
    print '+++ Velocity:\n',rnd(best_v,2)
    print '+++ Depth:\n',rnd(best_d,2)
    print '+++ VpVs:\n',rnd(best_r,2)

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

    init_plotting()

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

    ax1 = plt.subplot(221)

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

        if l == 'Min': ax1.plot(array(xs),-array(ys), linewidth=3, color='k', linestyle='--', label=l)
        elif l == 'Max': ax1.plot(array(xs),-array(ys), linewidth=3, color='k', linestyle='-.', label=l)
        else: ax1.plot(array(xs),-array(ys), linewidth=3, color=c, linestyle='-', label=l)
   
    ax1.set_xlabel('Velocity [km/s]')
    ax1.set_ylabel('Depth [km]')
    ax1.set_xlim(inp_dic['PLT_VP_RNG'])
    ax1.set_ylim(-array(inp_dic['PLT_DP_RNG'])[::-1])
    ax1.locator_params(axis = 'x', nbins = 6)
    ax1.locator_params(axis = 'y', nbins = 6)
    ax1.grid()
    ax1.legend(loc=3)

    ax2 = plt.subplot(224)

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

        if l == 'Min': ax2.plot(array(xs),-array(ys), linewidth=3, color='k', linestyle='--', label=l)
        elif l == 'Max': ax2.plot(array(xs),-array(ys), linewidth=3, color='k', linestyle='-.', label=l)
        else: ax2.plot(array(xs),-array(ys), linewidth=3, color=c, linestyle='-', label=l)

    ax2.set_xlabel('VpVs [km/s]')
    ax2.set_ylabel('Depth [km]')
    ax2.set_xlim(inp_dic['PLT_R_RNG'])
    ax2.set_ylim(-array(inp_dic['PLT_DP_RNG'])[::-1])
    ax2.locator_params(axis = 'x', nbins = 6)
    ax2.locator_params(axis = 'y', nbins = 6)
    ax2.grid()
    ax2.legend(loc=2)

    plt.tight_layout()
    plt.savefig(dbase_name+'.png')
    plt.close()


    #__________StdDev (V,D,R)
    
    models    = loadtxt('models.dat')
    model_std = std(models,axis=0)[0]
    model_men = mean(models,axis=0)[0]
    best      = array([best_v,best_d,best_r]).flatten()
    tot_ax    = ceil(models.shape[1]/4.)

    init_plotting()

    for par in range(models.shape[1]):

        ax = plt.subplot(tot_ax,4,par+1)
        #ax.set_title('par=%d'%par, fontsize=14)
        ax.text(0.97, 0.82, 'par=%d'%(par+1), fontsize=12,
                transform=ax.transAxes, bbox={'facecolor':'w', 'alpha':0.5, 'pad':2},
                horizontalalignment='right', verticalalignment='top')

        mu    = mean(models[:,par])
        sig   = std(models[:,par])
        x     = linspace(mu-3*sig,mu+3*sig, 100)
        a,b,c = plt.hist(models[:,par], 50, alpha=.6)
        plt.vlines(best[par],0,max(a),color='r',zorder=20,linewidth=2)
        plt.vlines(mu,0,max(a),color='y',zorder=20,linewidth=2)
        plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
        plt.locator_params(axis = 'x', nbins = 4)
        plt.locator_params(axis = 'y', nbins = 4)

    plt.tight_layout()
    plt.savefig('models_stat.png')
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
    plt.savefig(dbase_name+'_mut_r.png')
    plt.close()

best, best_rs, best_v, best_d, best_r = run_ga()
write_res(dbase_name, eval(resetDB), best_rs, best_v, best_d, best_r)
plot(best, best_rs, best_v, best_d, best_r)

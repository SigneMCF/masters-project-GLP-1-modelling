import argparse
import os 
import sys
sys.path.append("/mnt/code")
from cellsystem import CellSystem
from parameters_df import *

def main(args):
    myparam=Parameters()
    myparam.set("filename",args.filename)
    myparam.set("folder",args.folder)
    myparam.set("K_GR",args.K_GR)
    myparam.set("t_max",args.t_max)
    myparam.set("t_init",args.t_init)
    myparam.set("f_active",args.f_active)
    myparam.set("recycling_rate",args.recycling_rate)
    myparam.set("K_RL",args.K_RL)
    myparam.set("K_GRL",args.K_GRL)
    myparam.set("clearance_rate",args.clearance_rate)
    myparam.set("decay",args.decay)
    myparam.set("profile",args.profile)
    myparam.set("t_init",args.t_init)
    myparam.set("decay",args.decay)
    myparam.set("initdecay",args.initdecay)
    myparam.set("zeroint",args.zeroint)
    myparam.set("D",args.D)
    myparam.update()
    print("Running with")
    print(myparam.my_param)
    mysystem=CellSystem(parameters=myparam)
    mysystem.run()
    

if __name__ == "__main__":
    #assert numpyro.__version__.startswith("0.9.2")
    default=Parameters()
    
    parser = argparse.ArgumentParser(description="GPL1 diffusion")
    parser.add_argument("--filename",nargs="?",default="test.txt",type=str)
    parser.add_argument("--folder",nargs="?",default=default.get("folder"),type=str)
    parser.add_argument("--K_GR", nargs="?", default=default.get("K_GR"), type=float)
    parser.add_argument("--t_min", nargs="?", default=default.get("t_min"), type=int)
    parser.add_argument("--t_max", nargs="?", default=default.get("t_max"), type=int)
    parser.add_argument("--t_init", nargs="?", default=default.get("t_init"), type=int)
    parser.add_argument("--f_active", nargs="?", default=default.get("f_active"), type=float)
    parser.add_argument("--recycling_rate",nargs="?",default=default.get("recycling_rate"),type=float) #set to 10 times smaller than default=log(2)/3600. we change this back for testing says signe.
    parser.add_argument("--K_RL",nargs="?",default=default.get("K_RL"),type=float)
    parser.add_argument("--profile",nargs="?",default=default.get("profile"),type=str)# set profile default
    parser.add_argument("--K_GRL",nargs="?",default=default.get("K_GRL"),type=float)
    parser.add_argument("--clearance_rate",nargs="?",default=default.get("clearance_rate"),type=float)
    parser.add_argument("--decay",nargs="?",default=default.get("decay"),type=float)
    parser.add_argument("--initdecay",nargs="?",default=default.get("initdecay"),type=float)
    parser.add_argument("--zeroint",nargs="?",default=default.get("zeroint"),type=float)
    parser.add_argument("--D",nargs="?",default=default.get("D"),type=float)
    
    args = parser.parse_args()
    main(args)


# ---defaults---- 
#a_dict['K_RL']={'value': 10,'comment': 'Dissociation constant (nM) of poor, in-active binder: R+L <-> RL'}
#a_dict['K_GR']={'value': 1E+4,'comment': 'Dissociation constant (nM) of G+R <-> GR'}
#a_dict['K_GRL']={'value': 0.1,'comment': 'Dissociation constant (nM) of active binder: GR+L <-> GRL'}
#a_dict['f_active']={'value':0.1,'comment':'Fraction of maximally active receptors (GRL) compared to Rtot'}
#a_dict['r_tube']={'value': 100,'comment': 'Tube radius (nM) Intercellular distance/2'} #used to be called R
#a_dict['v']={'value': 0, 'comment': 'Advective flow velocity (nm/sec). Used to be 100 nm/sec.'}
#a_dict['D']={'value': 10*1E6,'comment': 'Ligand diffusion (nm^2/s)'}
#a_dict['decay']={'value': np.log(2)/(10*3600),'comment': 'Interstititum clearence rate (/sec). Default calculated from T_1/2=10 h'}
#a_dict['clearance_rate']={'value': np.log(2)/(1*3600),'comment': 'Receptor mediated clearance rate (/sec). Default calculated from T_1/2=1 h'}
#a_dict['recycling_rate']={'value': np.log(2)/(1*3600),'comment': 'Receptor recycling rate (/sec). Default calculated from T_1/2=1 h'}
#a_dict['k_d_fixed']={'value': np.log(2)/(2*3600),'comment': 'Fixed dissociation rate for active ligand-receptor binding (/sec). Default calculated from T_1/2= 2 h'}
#a_dict['V_reservoir']={'value': 1.0E+18,'comment': 'Reservoir/source volume (nm^3).'}
#a_dict['interstitium_vol_fraction']={'value': 0.01,'comment': 'Volume fraction of interstitium'}
#a_dict['r_min']={'value': 0.5E+6,'comment': 'Radius of source (nm)'}
#a_dict['r_max']={'value': 10.0E+6,'comment': 'Total simulation radius corresponding to distance from center of source (nm)'}
#a_dict['L0']={'value': 100,'comment': 'Initial source concentration (nM)'}
#a_dict['t_min']={'value': 1.0E-9,'comment': 'minimum time step in simulation (seconds)'}
#a_dict['t_max']={'value': 864000,'comment':'Total simulation time (seconds)'}
#a_dict['n_min']={'value': 100,'comment' : '# simulation steps at smallest time step'}
#a_dict['n_pr_10']={'value': 100,'comment': '# simulation steps pr decade (log-Time stepping)


#Vary decay, K_RL , K_GRL, clearance rate: 
# K_RL: x100, /100:
# K_GRL: x100, /100
#decay: x100, /100    (default = 1.925408834888737e-05
#clearance: x100, /100 (default = 1.925e-4


#running  simulations with:
#0: default
#1: K_RL=1000
#2: K_RL=0.1
#3: K_GRL=10
#4: K_GRL=0.01
#5: decay=1.925e-3
#6: decay=1.925e-7
#7: clearance=1.925e-2
#8: clearance=1.925e-6


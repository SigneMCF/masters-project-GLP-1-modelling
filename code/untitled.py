import argparse
from cellsystem import CellSystem
from parameters_df import *

def main(args):
    myparam=Parameters()
    myparam.set("filename",args.filename)
    myparam.set("save_to_folder",args.save_to_folder)
    myparam.set("K_GR",args.K_GR)
    myparam.set("t_max",args.t_max)
    myparam.set("f_active",args.f_active)
    myparam.set("recycling_rate",args.recycling_rate)
    myparam.set("K_RL",args.K_RL)
    myparam.set("K_GRL",args.K_GRL)
    myparam.set("clearance_rate",args.clearance_rate)
    myparam.set("decay",args.decay)
    myparam.update()
    mysystem=CellSystem(myparam)
    mysystem.run()
    
    
    
    #a_dict['decay']={'value': np.log(2)/(10*3600),'comment': 'Interstititum clearence rate (/sec). Default calculated from T_1/2=10 h'}
    #    a_dict['clearance_rate']={'value': np.log(2)/(1*3600),'comment': 'Receptor mediated clearance rate (/sec). Default calculated from T_1/2=1 h'}
    #    a_dict['recycling_rate']={'value': np.log(2)/(1*3600),'comment': 'Receptor recycling rate (/sec). Default calculated from T_1/2=1 h'}
  
   

#Vary decay, K_RL , K_GRL, clearance rate: 
# K_RL: x100, /100:
# K_GRL: x100, /100
#decay: x100, /100
#clearance: x100, /100
    
#myparam.set("K_RL",10000)
#myparam.set("K_GRL",100)
myparam.update()
myparam.my_param
    
    mysystem=CellSystem(parameters=myparam)
    


if __name__ == "__main__":
    #assert numpyro.__version__.startswith("0.9.2")
    default=Parameters()
    
    parser = argparse.ArgumentParser(description="GPL1 diffusion")
    parser.add_argument("-f","--filename",nargs="1",default="",type=str)
    parser.add_argument("--save_to_folder",nargs="1",default=default.get("save_to_folder",type=str)
    parser.add_argument("--K_GR", nargs="1", default=default.get("K_GR"), type=float)
    parser.add_argument("--t_max", nargs="1", default=default.get("K_GR"), type=int)
    parser.add_argument("--f_active", nargs="1", default=default.get("f_active"), type=float)
    parser.add_argument("--recycling_rate",nargs="1",default=default.get("recycling_rate"),type=float)
    parser.add_argument("--K_RL",nargs="1",nargs="1",default=default.get("K_RL"),type=float)
    parser.add_argument("--K_GRL",nargs="1",nargs="1",default=default.get("K_GRL"),type=float)
    parser.add_argument("--clearance_rate",nargs="1",default=default.get("clearance_rate"),type=float)
    parser.add_argument("--decay",nargs="1",default=default.get("decay"),type=float)
    
    args = parser.parse_args()
    main(args)

    #numpyro.set_platform(args.device)
    #numpyro.set_host_device_count(args.num_chains)

    #main(args)
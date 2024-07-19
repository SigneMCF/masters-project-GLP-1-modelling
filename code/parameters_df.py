import numpy as np
import pandas as pd


#myparam.set("recycling_rate",1.925409e-04/10)

class Parameters:
    
    Mol_times_nm3=0.6022
    
    def __init__(self):
        a_dict=dict()
        a_dict['r_cell']={'value': 10000,'comment': 'Cell radius in nm'} 
        a_dict['b']={'value': 10,'comment': 'Surface thickness in nm'}
        a_dict['K_RL']={'value': 10,'comment': 'Dissociation constant (nM) of poor, in-active binder: R+L <-> RL'}
        a_dict['K_GR']={'value': 1E+4,'comment': 'Dissociation constant (nM) of G+R <-> GR'}
        a_dict['K_GRL']={'value': 0.1,'comment': 'Dissociation constant (nM) of active binder: GR+L <-> GRL'}
        a_dict['f_active']={'value':0.1,'comment':'Fraction of maximally active receptors (GRL) compared to Rtot'}
        a_dict['r_tube']={'value': 100,'comment': 'Tube radius (nm) Intercellular distance/2'} #used to be called R
        a_dict['v']={'value': 0, 'comment': 'Advective flow velocity (nm/sec). Used to be 100 nm/sec.'}
        a_dict['D']={'value': 10*1E6,'comment': 'Ligand diffusion (nm^2/s)'}
        a_dict['decay']={'value': np.log(2)/(10*3600),'comment': 'Interstititum clearence rate (/sec). Default calculated from T_1/2=10 h'}
        a_dict['initdecay']={'value': 0.0014440566261665528,'comment': 'initial decay for diffusion, unit/ sec'}
        a_dict['clearance_rate']={'value': np.log(2)/(1*3600),'comment': 'Receptor mediated clearance rate (/sec). Default calculated from T_1/2=1 h'}
        a_dict['recycling_rate']={'value': np.log(2)/(1*3600),'comment': 'Receptor recycling rate (/sec). Default calculated from T_1/2=1 h'}
        a_dict['k_d_fixed']={'value': np.log(2)/(2*3600),'comment': 'Fixed dissociation rate for active ligand-receptor binding (/sec). Default calculated from T_1/2= 2 h'}
        a_dict['V_reservoir']={'value': 1.0E+18,'comment': 'Reservoir/source volume (nm^3).'}
        a_dict['interstitium_vol_fraction']={'value': 0.01,'comment': 'Volume fraction of interstitium, used to calculate volume of receptors'}
        a_dict['r_min']={'value': 0.5E+6,'comment': 'Radius of source (nm)'}
        a_dict['r_max']={'value': 10.0E+6,'comment': 'Total simulation radius corresponding to distance from center of source (nm)'}
        a_dict['L0']={'value': 100,'comment': 'Initial source concentration (nM)'}
        a_dict['zeroint']={'value': 1,'comment': 'shuts of internalisation'}
        a_dict['t_min']={'value': 1.0E-9,'comment': 'minimum time step in simulation (seconds)'}
        a_dict['t_max']={'value': 864000,'comment':'Total simulation time (seconds)'}
        a_dict['t_init']={'value': 3600,'comment':'Time passed at the start of simulation(seconds)'}
        a_dict['n_min']={'value': 100,'comment' : '# simulation steps at smallest time step'}
        a_dict['n_pr_10']={'value': 100,'comment': '# simulation steps pr decade (log-Time stepping)'}
        a_dict['folder']={'value': 'Results','comment': 'Pass a folder to save results.'}# altered for sensitivity analysis
        a_dict['filename']={'value': 't_Lall.txt','comment': 'Filename for results (if a folder has been given)'}
        a_dict['n_t_saved']={'value' : 100,'comment': '# of time points saved'}
        a_dict['profile']={'value' : "diffusion",'comment': 'baseline profile'}
    
        self.my_dict=a_dict
        keys=a_dict.keys()
        values=[x['value'] for x in a_dict.values()]
        comments=[x['comment'] for x in a_dict.values()]
        self.my_param=pd.DataFrame({'value': values,'comment': comments},index=keys)
        
        #calculated properties
        self.calc_dict=dict()
        self.calc_dict['Rtot']={'value': None,'comment': 'Total density of receptors (nM)'}
        self.calc_dict['k_d_poor']={'value': None,'comment': 'Dissociation rate of receptor-inactive ligand binding (/sec)'}
        self.calc_dict['k_a']={'value': None,'comment': 'Association rate for receptor-ligand binding (/(sec*nM)). Assumed the same for both active and inactive receptor'}
        self.calc_dict['K_RL_clearance']={'value': None,'comment': 'Effective dissocation constant for receptor-inactive ligand binding, R+L <-> RL (nM), shifted due to clearance.'}
        self.calc_dict['Gtot']={'value': None,'comment': 'Total concentration of G (nM)'}
    
        self.calc_dict['cooperativity']={'value': None,'comment': 'Cooperativity of GR-> GRL (RL -> GRL) compared to R->RL (R -> GR)'}
 
        self.calc_dict['epsilon']={'value': None,'comment': 'Reaction volume at cell surface relative to interstitium volume'}
        self.update()
    
    def get(self,key):
        if key in self.calc_dict.keys():
            val=self.calc_dict[key]['value']
        else:
            val=self.my_param.loc[key]['value']
        try:
            fval=float(val)
        except:
            fval=val
        return fval
    
    def get_comment(self,key):
        if key in self.calc_dict.keys():
            return self.calc_dict[key]['comment']      
        return self.my_param.loc[key]['comment']
    
    def update(self):
        Rtot=50000/(4*np.pi*np.power(self.get('r_cell'),2)*self.get('b'))*self.Mol_times_nm3*1E9 #density of receptors (nM)
        k_d_poor=self.get('K_RL')/self.get('K_GRL')*self.get('k_d_fixed')
        k_a=self.get('k_d_fixed')/self.get('K_GRL')
        K_RL_clearance=(k_d_poor+self.get('clearance_rate'))/k_a
        self.set('Rtot',Rtot)
        self.set('kd_poor',k_d_poor)
        self.set('k_a',k_a)
        self.set('K_RL_clearance',K_RL_clearance)
      
        alpha=self.get('K_RL')/self.get('K_GRL')
        self.set('cooperativity',alpha)
        self.set('epsilon',2*self.get('b')/self.get('r_tube') )
        
        #formular for K from binary reaction,    
        f_active=self.get('f_active')
        Gtot=f_active*Rtot
        self.set('Gtot',Gtot)
        #K_GR=Rtot*f_active+Gtot/f_active-(Rtot+Gtot)
        #self.set('K_GR',K_GR)
         
    def set(self,key,value):
        if key in self.calc_dict:
            self.calc_dict[key]['value']=value
        elif key in set(self.my_param.index):
            self.my_param.loc[key,'value']=value
            
    def set_parameters(self,param):
        self.my_param=param
        self.update()
    
    def get_parameters(self):
        return self.my_param
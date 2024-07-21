# -*- coding: utf-8 -*-

import numpy as np
import matplotlib.pyplot as plt
import scipy
from scipy.special import jv 
from scipy.sparse import diags
from numpy.linalg import inv
from scipy.integrate import odeint
from scipy.interpolate import InterpolatedUnivariateSpline
from scipy.integrate import quad
from scipy import special
import pandas as pd
from parameters_df import Parameters
from cellmodel import CellModel
import os

debug=0
# old diffusion
# diffusion integrator for initial state. this one has a faulty decay term.
#def diffusion2(x,d,t,c0,x0,loss):
#    return(c0*np.exp(-loss*t)*(1-special.erf((x-x0)/(np.sqrt(4*d*t)))))
# current diffusion
def diffusion2(x,d,t,c0,x0,loss):
    u0=(x-x0)/np.sqrt(2*d*t)
    res=np.ones(len(x))

    for i in range(len(x)):
        def inting(u):
             return(np.sqrt(2/np.pi)*np.exp(-(u**2)/2 -(u**-2)*(1/2)*((x[i]-x0)/(np.sqrt(d/loss)))**2))
        res[i]=quad(inting,u0[i],15)[0]*c0
    return res

# integral of diffusion, for if we wanted entire segments and not just points. without loss.
def integrator(c0,d,t,xvec):
    totallength=0
    cells=np.zeros(len(xvec))
    for i in range(len(xvec)):
        cells[i]=quad(diffusion2,totallength,totallength+xvec[i],args=(d,t,c0))[0]
        totallength+=xvec[i]
    return cells

# this is where log timing start, we can fiddle with this for linear time.

def LogSpacing(first,last,step_pr_decade):
    if debug==1:
        print("running logspacing")
    ndecades=round(np.log(last/first)/np.log(10))
    ns=int(round(ndecades*step_pr_decade+1))
    return np.exp(np.linspace(np.log(first),np.log(last),ns)) #round(ndecades*step_pr_decade+1)))

def Const_LogSpacing(first,last,step_pr_decade,repeats):
    if debug==1:
        print("running constant logspacing")
    equil=np.linspace(first,first*repeats,int(repeats))
    new_first=equil[-1]+first
    if (new_first>last):
        return equil
    return np.concatenate((equil,LogSpacing(new_first,last,step_pr_decade))) #used for time initializations
    
# this is for the old variant with advection.
def half_length(p):
    return np.log(2)*2*p.D/(np.sqrt(p.v*p.v+4*p.D*p.decay)-p.v)
    
# these are our gradients and diffs. Note that some of these tech
def Divergence(r,v):
    v_sp=InterpolatedUnivariateSpline(r, v, k=2)
    d_v=v_sp.derivative()
    return d_v(r)+2/r*v

def Derivative(r,v,n=1):
    v_sp=InterpolatedUnivariateSpline(r, v, k=2)
    d_v=v_sp.derivative(n)
    return d_v(r)

def Integrate(r,v):
    r2=r*r
    v_sp=InterpolatedUnivariateSpline(r, r2*v, k=2)
    return 4*np.pi*v_sp.integral(r[0],r[-1])


def Laplace(r,v):
    v_sp=InterpolatedUnivariateSpline(r, v, k=2)
    d_v=v_sp.derivative()
    d_v2=v_sp.derivative(2)
    return 2/r*d_v(r)+d_v2(r)
# finds indices
def find_indices(avec,values):
    indices=np.zeros([len(values)],dtype=int)
    for i in range(0,len(values)):
        indices[i]=np.argmin(abs(avec-values[i]))
    return indices

# defines the states etc.
class CellSystem:
    
    def __init__(self,parameters):
        self.P=parameters
        self.cellmodel=CellModel()
        self.cellmodel.set(parameters)
        self.cellmodel.set_geometry('bulk_and_surface',parameters.get('epsilon'))
        
        #FEM:
        
        self.N_t=None #number of discretizations in time
        self.N_r=200  #number of steps in radial direction
        
        #Internal registration 
        self._current_t_index=None
        self._file_obj=None
    
        #Solution vectors
        self.T=None
        self.r=None     
        self.a=None
        self.ltot=None
        self.Ltot_inflow=None
        
        self.R=None
        self.RL=None
        self.GR=None
        self.GRL=None
        self.L0=None #start concentration in nM at r=r_min
        self.Ltot=None
                 
        #left and right boundary , as defined by the normal direction (-1 and +1)
        self.boundary=None
        self.set_boundary()
 
    def _open_file(self):
        folder=self.P.get('folder')
        filename=self.P.get('filename')
        if len(folder)>1 and os.path.exists(folder) and os.path.isdir(folder) and len(filename)>1:
            fullname=os.path.join(folder,filename)
            if os.path.exists(fullname):
                os.remove(fullname)
            self._file_obj=open(fullname,'a')
        else:
            self._file_obj=None
# to write our results    
    def _write_to_file(self,t,Lall):
        if not self._file_obj is None:
            Lall_txt='\t'.join(np.array(Lall,'str'))
            line=str(t)+'\t'+Lall_txt+'\n'
            self._file_obj.write(line)  
                
    def _close_file(self):
        if self._file_obj is not None:
            self._file_obj.close()
    # to read our results, remember we only save the Lall, and calculate everything from that
    def read_file(self,file,solve_for_cells=True):
        with open(file,"r") as f:
            content=f.readlines()
        content=np.array([x.split('\t') for x in content],float)
        s=content.shape
        self.N_r=s[1]-2 #first index is time, next is Ltot_inflow
        T=content[:,0] 
        
        self._set_FEM(T)     
        self._allocate()
        
        self.Ltot_inflow=content[:,1]
        self.Ltot=content[:,2:]
        
        self._initialize_cells(from_tot=True)
        
        if solve_for_cells:
            self._solve_for_cells()

     # get our simulation distances and makes segments   
    def get_r(self):
        return np.linspace(self.P.get('r_min'),self.P.get('r_max'),self.N_r)
    # gets all our parameters and the geometry
    def set_parameters(self,parameters):
        self.P=parameters
        self.cellmodel.set(parameters)
        self.cellmodel.set_geometry('bulk_and_surface',parameters.get('epsilon'))
 

    # here be the boundary conditions yarrr,  for future work this is where the decaying source could be implemented.

    def set_boundary(self,type="fixed"):
        if type=="reservoir" or type=="fixed":
            self.boundary=[{'type': type,'value': self.L0,'normal': -1},{'type': 'value','value':0,'normal':1}]
        elif type=="filter":
            self.boundary=[{'type':'filter','value':0,'normal': -1},{'type': 'filter','value': 0,'normal': 1}]
        elif type=="reflective":
            self.boundary=[{'type': 'gradient','value':0,'normal': -1},{'type': 'gradient','value': 0,'normal': 1}]
        else:
            print("Unrecognized boundary type")
            
    # sets some time and space things
    def _set_FEM(self,T,r_vector=None):
        self.T=T
        self.N_t=len(T)
        if r_vector is None: 
            self.r=self.get_r()
        else:
            self.r=r_vector
            self.P.set('r_min',r_vector[0])
            self.P.set('r_max',r_vector[-1])
            self.N_r=len(r_vector)
    
    def _allocate(self):
        # allocate, including some old variants.

        self.L=np.zeros([self.N_t,self.N_r])
        #self.L[0,:]+=initconc
        #self.L=np.zeros([self.N_t,self.N_r]) #concentration of ligands in interstitium = free ligand concentration at cell surface
        #self.Ltot=np.zeros([self.N_t,self.N_r]) #total concentration of ligands = L+RL+GRL
        self.Ltot_inflow=np.zeros([self.N_t])
        #self.rho0=np.zeros([self.N_t,self.N_z])
        #self.rho0[0,:]=self.rho
        self.R=np.zeros([self.N_t,self.N_r])
        self.RL=np.zeros([self.N_t,self.N_r])
        self.GR=np.zeros([self.N_t,self.N_r])
        self.GRL=np.zeros([self.N_t,self.N_r])
        self.Rint=np.zeros([self.N_t,self.N_r])
        self.Ltot=np.zeros([self.N_t,self.N_r])     
        self._last_solutions=np.zeros([self.N_r,len(self.cellmodel.molecule_enum)])
        

        
        self.init_solutions=np.zeros([self.N_r,len(self.cellmodel.molecule_enum)])

        
       
        
        self.integrated_GRL=np.zeros([self.N_t]) 
    # initial condition profiles, most of these are old and not used.
    def _set_profile(self,profile):
        if 'L0' in profile:
            self.L0=profile['L0'] #start concentration in nM at r=r_min
        else:
            self.L0=self.P.get('L0')
        #print(self.l_i0)
 
        if profile['type']=='flat':
            self.L[0,:]=self.L0  
        #put densities linear decaying from z=0 to z=1  
        elif profile['type']=='linear':       
            for zi in range(0,self.N_r):
                self.L[0,zi]=(1-zi/(self.N_r-1))*self.L0  #put densities linear decaying from z=0 to z=1
                                
        elif profile['type']=='dirac':
            self.L[0,0]=self.L0
        elif profile['type']=='soft':
            x0=round(profile['length']*self.N_r)  #length is given as a fraction of the total length. x0 are corresponding # indices 
            #print('Length={0} l_i0={1}'.format(x0,self.l_i0))
            M=np.array([[x0,x0*x0/2],[x0*x0/2,x0*x0*x0/6]])
            invM=inv(M)
            a=np.dot(invM,np.array([0,-self.L0]))
            indices=np.arange(0,x0)
            indices2=indices*indices
            indices3=indices2*indices
            l=self.L0+a[0]/2*indices2+a[1]/6*indices3
            #print(l)
            self.L[0,indices]=l
        elif profile['type']=='last':
            self.L[0,:]=self.last_L
        elif profile['type']=='diffusion':
            #xvec=np.ones(self.N_r)*round(self.P.get("r_max")/self.N_r)
            #initconc=integrator(self.P.get("L0"),self.P.get("D"),self.P.get("t_init"),xvec)
            initconc=diffusion2(self.r,self.P.get("D"),self.P.get("t_init"),self.L0,self.r[0],self.P.get("initdecay"))
            self.L[0,:]=initconc
        else:
            print("Cannot recognize "+profile['type']+" as a profile type: Options are flat,linear,soft, dirac or last")
            return False                                
        
    # intialization
        
    def _initialize(self,T,r_vector,profile): #T is end time , r_vector is vector of radii for simulation
 
        self._set_FEM(T,r_vector)
        self._allocate()
        
        self._current_t_index=0
        self._set_profile(profile)                       
        
        self._initialize_cells(from_tot=False)
        
        #self.l_c[0,:]=self.calc_l_c(self.ltot[0,:])
        self._open_file()
        return True
    
    
    
    

        


# initilization for each segment
    def _initialize_cells(self,from_tot=False):
        last_solution=None

        
        
        for i in range(self.N_r-1,-1,-1):
            #print("Doing initialization for i={}".format(i))
            if from_tot:
                self.cellmodel.Ltot=self.Ltot[0,i]
                self.cellmodel.L=None
            else:
                self.cellmodel.L=self.L[0,i]
                self.cellmodel.Ltot=None
            #print("Init cells, doing i={}".format(i))
            last_solution=self.cellmodel.solve(last_solution)
            self.init_solutions[i,:]=np.ndarray.copy(last_solution)
            self.R[0,i]=self.cellmodel.get('R')

            if from_tot:
                self.L[0,i]=self.cellmodel.get('L')
            else:
                self.Ltot[0,i]=self.cellmodel.get('Ltot')
           
            self.RL[0,i]=self.cellmodel.get('RL')
            self.GR[0,i]=self.cellmodel.get('GR')
            self.GRL[0,i]=self.cellmodel.get('GRL')
            self.Rint[0,i]=self.cellmodel.get('Rint')
        
        self._last_solutions=np.ndarray.copy(self.init_solutions)
        self.cellmodel.L=None 
    # solving the segments
    def _solve_for_cells(self,at_indices=None):
        if debug==1:
            print("running solve for cells")
        if at_indices is None:
            at_indices=range(0,self.N_t)
        #epsilon is reaction_volume relative to interstitium volume
        eps=self.P.get('epsilon')
        receptor_vol_frac=eps*self.P.get('interstitium_vol_fraction')
        
        #fill out values:
        #print("Done, filling out cell concentrations...")
        self._last_solutions=np.ndarray.copy(self.init_solutions)
        #print("shape of Lall={}".format(np.shape(Lall)))
        print("# time indices ={}".format(len(at_indices)))
        for ti in at_indices:
            if (ti % 100) ==0:
                print("doing ti={}".format(ti))
                
            for zi in range(self.N_r-1,-1,-1):
                self.cellmodel.Ltot=max(0,self.Ltot[ti,zi])
                self._last_solutions[zi,:]=np.ndarray.copy(self.cellmodel.solve(self._last_solutions[zi,:]))
                self.R[ti,zi]=self._last_solutions[zi,self.cellmodel.iof('R')]
                self.L[ti,zi]=self._last_solutions[zi,self.cellmodel.iof('L')]
                self.RL[ti,zi]=self._last_solutions[zi,self.cellmodel.iof('RL')]
                self.GR[ti,zi]=self._last_solutions[zi,self.cellmodel.iof('GR')]
                self.GRL[ti,zi]=self._last_solutions[zi,self.cellmodel.iof('GRL')]
                
            #self.R[ti,:]=np.ndarray.copy(self._last_solutions[:,self.cellmodel.iof('R')])
            #self.L[ti,:]=np.ndarray.copy(self._last_solutions[:,self.cellmodel.iof('L')])
            #self.RL[ti,:]=np.ndarray.copy(self._last_solutions[:,self.cellmodel.iof('RL')])
            #self.GR[ti,:]=np.ndarray.copy(self._last_solutions[:,self.cellmodel.iof('GR')])
            #self.GRL[ti,:]=np.ndarray.copy(self._last_solutions[:,self.cellmodel.iof('GRL')])
            self.integrated_GRL[ti]=receptor_vol_frac*1E-24*Integrate(self.r,self.GRL[ti,:]) #1E-24 to get from (nm^3)-> liter
        self._solved_for_cells=True
    


  # this is where diffusion and stuff happens, aka the actual differential part.
    
    def dLall_dt(self,Lall,t):
        self._number_of_calls=self._number_of_calls+1
        
        #L=np.zeros(self.N_r)
        #RL=np.zeros(self.N_r)
        Lall[Lall<0]=0
        Ltot_inflow,Ltot=Lall[0],Lall[1:]
        at_which_index=np.argmax(self.T>t)
        if at_which_index>self._current_t_index+1:
            self._current_t_index=at_which_index-1
            print(f"Finished time {t} corresponding to index {self._current_t_index}")
            self._write_to_file(t,Lall)
            #write to result file
        #print('#calls={0}, t={1}'.format(self._number_of_calls,t))
        #if sum(Lall<0)>0:
        #    print("Warning - negative Lall values")
        #    Lall[Lall<0]=0
            
        Ltot_inflow,Ltot=Lall[0],Lall[1:]
        for zi in range(0,self.N_r):
            self.cellmodel.Ltot=Ltot[zi]
            self._last_solutions[zi,:]=self.cellmodel.solve(self._last_solutions[zi,:])
            
        
        L=np.ndarray.copy(self._last_solutions[:,self.cellmodel.iof('L')]) #get free ligand concentration from Ltot
        RL=np.ndarray.copy(self._last_solutions[:,self.cellmodel.iof('RL')])
       
    # what are the div and lap functions here
        div_L=Divergence(self.r,L)
        lap_L=Laplace(self.r,L)
        d_L=Derivative(self.r,L)
        
        eps=self.P.get('epsilon')  
        D=self.P.get('D')
        v=self.P.get('v')
        decay=self.P.get('decay')
        clearance=self.P.get('clearance_rate')
    
        dLtot_dt=-v*div_L+D*lap_L-decay*L-eps*clearance*RL
    
        j0=v*L[0]-D*d_L[0]     
        A=np.pi*np.power(self.P.get('r_tube'),2)
        if self.boundary[0]['type']=='reservoir': 
            dLtot_dt[0]=-A/self.P.get('V_reservoir')*j0
        elif self.boundary[0]['type']=='fixed':
            dLtot_dt[0]=0
        dLtot_inflow_dt=A*j0*1E-24 #nM*(nm)^3/s -> nmol/s
        #print("dLtot_inflow_dt={}".format(dLtot_inflow_dt))
        self._last_Lall=Lall
        return np.concatenate(([dLtot_inflow_dt],dLtot_dt))
    
    #This is the main function to call:
    
    def solve(self,T,r_vector=None,profile=None,solve_for_cells=True):
        
        self.P.update()
        self._solved_for_cells=False
        ok=self._initialize(T=T,r_vector=r_vector,profile=profile)
        if not ok:
            return
      
        Ltot_inflow=self.Ltot_inflow[0]
        Ltot=self.Ltot[0,:]
        Lall=np.concatenate(([Ltot_inflow],Ltot))
        #print('l0={} '.format(l0))
        #print(T)
        self._number_of_calls=0
        Lall=odeint(self.dLall_dt,Lall,t=self.T)
        self._close_file()
        
        #epsilon is reaction_volume relative to interstitium volume
        eps=self.P.get('epsilon')
        receptor_vol_frac=eps*self.P.get('interstitium_vol_fraction')
        
        #fill out values:
        #print("# time indices ={}".format(self.N_t))
        for ti in range(0,self.N_t):
            self.Ltot[ti,:]=Lall[ti,1:]
            self.Ltot_inflow[ti]=Lall[ti,0]
            
        if solve_for_cells:
            self._solve_for_cells()
            
        self.last_L=np.ndarray.copy(self.L[-1,:])
    # this is to run things
    def run(self,T=None): #init_type='diffusion'
        if debug==1:
            print("running run,heh")
        if T is None:
            t_min=self.P.get('t_min')
            t_max=self.P.get('t_max')
            n_min=self.P.get('n_min')
            n_pr_10=self.P.get('n_pr_10')
            T=Const_LogSpacing(t_min,t_max,n_pr_10,n_min)
        self.set_boundary("fixed")
        self.N_r=600
        profile={'type': self.P.get("profile")} #default used to be 'soft'
        profile['length']=0.025
        self.solve(T=T,profile=profile,solve_for_cells=False)
    
    # indicing
    def find_indices(self,hours=None):
        if hours is None:
            return None
        tr=find_indices(self.T,hours*3600)
        return tr
    # some built in plots, these are used primarily to get an understanding of the long time model equilibria     
    def plot(self,indices=None):
        if debug==1:
            print("running plot")
            
        myr=self.r/1000 #in micron
        print("length of simulation region in micron is {}".format(len(myr)))
        plt.figure(figsize=(12,18)) #width,height in inches
        lT=len(self.T)
        #print('Length of t {}'.format(self.T))
        if indices is None:
            tr=np.arange(0,lT,1)
        else:
            tr=indices
        #print(tr)
        if not lT-1 in tr:
            tr=np.concatenate((tr,np.array([lT-1])))
        #print(tr)
        if not self._solved_for_cells:
            self._solve_for_cells(tr)
        
        T_in_h=np.array(self.T)/(60*60)
        #print(T_in_h)
        for plot_nr in range(0,6):
            plt.subplot(3,2,plot_nr+1)
            if plot_nr<5:
                if plot_nr==0:
                    #print("li shape {}".format(self.l_i.shape))
                    q=self.L
                    tit=r'$L$'
                if plot_nr==1:
                    #print("rho1 shape {}".format(self.rho1.shape))
                    q=self.RL
                    tit=r'$RL$'
                if plot_nr==2:
                    q=self.GR
                    tit=r'$GR$'
                if plot_nr==3:
                    q=self.GRL
                    tit=r'$GRL$'
                if plot_nr==4:
                    #q=self.Ltot
                    #tit=r'$Ltot$'
                    q=self.R+self.RL+self.GRL+self.GR
                    tit='$R_t$'
                    
                #print(q.shape)
                for t in tr:
                    plt.plot(myr,q[t,:],label="t={:.2E} h".format(T_in_h[t]))
            
                plt.title(tit)
                plt.xlabel(r'r ($\mu$ m)')
                plt.ylabel(r'conc. (nM)')
                plt.legend(loc='upper right')
                
            #elif plot_nr==2:
                #tit=r'$L_t$'
                #plt.plot(T_in_h,self.Ltot_inflow)
                #plt.title(tit)
                #plt.xlabel(r't (h)')
                #plt.xscale('log')
                #plt.ylabel(r'n (nmol)')
            else:
                plt.figure(figsize=(10,6))
                plt.title('Total GRL')
                plt.plot(T_in_h[tr],self.integrated_GRL[tr])
                plt.xlabel(r't (h)')
                plt.xscale('log')
                plt.ylabel(r'$N(GRL)$ (nmol)')
        
        plt.show()
        
# the last parts are mostly outdated and were used for bugfixing and old runs .    

   # def get_last(self):
   #     val_dict=dict()
   #     unit_dict=dict()   
   #     tot_bound=self.integrated_rho1[-1]
   #     val_dict['Tot_bound']=tot_bound
    #    unit_dict['Tot_bound']='nmol'
        
    #    normed_rho1=self.rho1[-1,:]/(self.P.get('rho')* self.P.get('fraction')) 
    #    fracs=np.array([0.1,0.25,0.5])
     #   myr=self.r/1000
       
    #    for fr in fracs:
    #        i=np.argmin(np.abs(normed_rho1-fr))
    #        key="r_{}".format(fr)
    #        val_dict[key]=myr[i]
    #        unit_dict[key]='Micron'
    #    
    #    return val_dict,unit_dict
   
    #def print_last(self):
    #    val_dict,unit_dict=self.get_last()
     #   print('At time {}:'.format(self.T[-1]/3600))
     #   print('-------------------:')
     #   for key in val_dict.keys():
      #      print('{0}: {1} = {2} ({3})'.format(key,val_dict[key],unit_dict[key]))
    
def generate_scan(**kwargs):
    if debug==1:
        print("running scan")
    #kwargs is a dictionary. For each key is a vector of values to be tested. We run through all combinations
    result_df=None
    n_keys=len(kwargs)
    n_elements=np.zeros(n_keys,dtype=int)
    for i,key in enumerate(kwargs):
        n_elements[i]=len(kwargs[key])
    
    size_elements=np.zeros(n_keys,dtype=int)
    total=1        
    for i,key in enumerate(kwargs):
        size_elements[-(i+1)]=total
        total=total*n_elements[-(i+1)]
    
    #print(size_elements)
    ret_df=pd.DataFrame()
    for key in kwargs.keys():
        ret_df[key]=np.zeros(total,dtype=type(kwargs[key]))
    for all_combination in range(0,total):
        indices=np.zeros(n_keys,dtype=int)
        all_val=all_combination
        for i in range(0,n_keys):
            indices[i]=all_val//size_elements[i]
            all_val=all_val % size_elements[i]
      
        #print(indices)   
        for i,key in enumerate(kwargs):
            ret_df.loc[all_combination,key]=kwargs[key][indices[i]]
    
    return ret_df
            

#parameters is a dataframe with a row for each parameter combination and a column for each parameter to be changed Parameters can be generated with generate_scan
def generate_solutions(parameters,t_max=3600*24*20,n_step_pr_10=1000,file=None):
    if debug==1:
        print("running generate solutions")
    #print(parameters)
    keys=parameters.columns
    results=None
    t_min=1E-9
    n_min=2000
    #t_max=3600*24*20
    #t_max=1E-7
    #n_step_pr_10=1000
    myT=Const_LogSpacing(t_min,t_max,n_step_pr_10,n_min)
    
    for i in range(0,len(parameters)):
        param_i=Parameters()
        for k in keys:
            param_i.set(k,parameters.loc[i,k])
        Solver=Glp_R_pde(parameters=param_i)
        Solver.set_boundary("fixed")
        Solver.model="simple"
        Solver.N_r=600
        Solver.solve(T=myT,profile={'type': 'soft','length':0.025})
        
        #indices=[int(x) for x in np.linspace(n_min,len(myT)-1,10)]
        #Solver.plot(indices=indices)
        #Solver.print_last()
        last_vals,units= Solver.get_last()
        #print(last_vals)
        if results is None:
            results=parameters
            for k in last_vals.keys():
                results[k]=np.zeros(len(parameters))
                results['{}_unit'.format(k)]=units[k]
        for k in last_vals.keys():
            results.loc[i,k]=last_vals[k]
        print('At step {}:'.format(i))
        print(results.loc[i])
        print('-------------------')
        if not file is None:
            results.to_csv(file,index=False)
    
    return results
            
        


            
        

        
        
    

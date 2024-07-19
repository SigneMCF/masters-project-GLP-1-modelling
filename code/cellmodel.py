import numpy as np
from scipy.optimize import fsolve

class CellModel:
    
    def __init__(self):
        self.molecule_enum={'L':0,'R':1,'G':2,'RL':3,'GR':4,'GRL':5}
        self.Ltot=1 #set Ltot to a given value and L to None, if constraint is given in terms of Ltot 
        self.L=None #set L to a given value and Ltot to None, if constraint is given in terms of L 
        self.Rtot=1
        self.Gtot=1
        self.K_RL=1
        self.K_GR=1
        self.alpha=1
        self.K_RL_clearance=1
        self.kin_over_kout=0
        self._result=None
        #self.rho_min=1
        self.set_geometry()
    
    def iof(self,name): #index of molecule
        return self.molecule_enum[name]
    
    def set_geometry(self,type='surface_only',epsilon=0.1):
        if type=='surface_only': 
            self.bulk=0;self.epsilon=1
        elif type=='bulk_and_surface' or type=='surface_and_bulk':
            self.bulk=1;self.epsilon=epsilon
        elif type=='bulk_only':
            self.bulk=1;self.epsilon=0    
    
    def set(self,parameters):
        self.K_RL=parameters.get('K_RL')
        self.K_GR=parameters.get('K_GR')
        self.alpha=parameters.get('cooperativity')
        self.K_RL_clearance=parameters.get('K_RL_clearance')
        kout=parameters.get('recycling_rate')
        kin=parameters.get('clearance_rate')
        self.zeroint=parameters.get('zeroint')
        if kout>0 and kin>0:
            self.kin_over_kout=kin/kout
        else:
            self.kin_over_kout=0
            
        self.Rtot=parameters.get('Rtot')
        self.Gtot=parameters.get('Gtot')
        K_min=np.min([self.K_RL,self.K_GR])
        #self.rho_min=np.min([K_min,K_min/self.alpha,self.G_tot,self.R_tot])
    
    def get_constraints(self,x):
        l=self.molecule_enum['L']
        r=self.molecule_enum['R']
        g=self.molecule_enum['G']
        rl=self.molecule_enum['RL']
        gr=self.molecule_enum['GR']
        grl=self.molecule_enum['GRL']
        K_GRL=self.K_RL/self.alpha
        #rin=6 rin=kin/kout*x[rl]
        f=np.zeros(6)
        rin=x[rl]*self.kin_over_kout*self.zeroint
        if self.L is None:
            f[l]=self.bulk*x[l]+self.epsilon*(x[l]+x[rl]+x[grl])-self.Ltot
        else:
            f[l]=x[l]-self.L
        f[r]=x[r]+x[rl]+x[gr]+x[grl]+rin-self.Rtot
        f[g]=x[g]+x[gr]+x[grl]-self.Gtot
        f[rl]=x[rl]- x[r]*x[l]/self.K_RL_clearance
        f[gr]=x[gr]-x[r]*x[g]/self.K_GR
        f[grl]=x[grl]-x[gr]*x[l]/K_GRL
        return f

    def solve(self,x_guess=None):
        l=self.molecule_enum['L']
        r=self.molecule_enum['R']
        g=self.molecule_enum['G']
        rl=self.molecule_enum['RL']
        gr=self.molecule_enum['GR']
        grl=self.molecule_enum['GRL']
        if x_guess is None:
            x_guess=np.zeros(len(self.molecule_enum))
        if self.L is not None:
            x_guess[self.iof('L')]=self.L
        self._result=fsolve(self.get_constraints,x_guess)
        self._result[self._result<0]=0
        return self._result

    def _get(self,name):
        return self._result[self.molecule_enum[name]]
    
    def get(self,name):
        if name=='Rint':
            return self._get('RL')*self.kin_over_kout
        elif name=='Ltot':
            return self.bulk*self._get('L')+self.epsilon*(self._get('L')+self._get('RL')+self._get('GRL'))
        elif name=='Gtot':
            return self._get('G')+self._get('GR')+self._get('GRL')
        elif name=='Rtot':
             return self._get('R')+self._get('RL')+self._get('GR')+self._get('GRL')+self._get('RL')*self.kin_over_kout
    
        return self._get(name)
    
    

    
    
    
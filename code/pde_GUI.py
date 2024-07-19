
import ipywidgets as widgets
import qgrid
from IPython.display import display, clear_output
import pandas as pd
from parameters_df import *
from cellsystem import *

class pde_GUI:
    
    def __init__(self):
        self.param=Parameters()
        self.mySystem=None
        self._param_grid = qgrid.QgridWidget(df=self.param.get_parameters())
        self._output=widgets.Output()
        self._output_box = widgets.Box(children=[self._output],
                                   layout={'height': '180px', 'border': '1px solid black', 'overflow': 'scroll'})
        
        self._run_button = widgets.Button(description='Run')
        self.my_ui=widgets.VBox([self._param_grid,self._run_button,self._output_box])
        self._set_callback()

 
        #with self.param_out:
        #    clear_output()
        #    self.param_grid.df=param_df 
        
        
        self._param_grid.on('cell_edited',self.cell_edit)
  
    def cell_edit(self,change, param_grid):
        param_df=param_grid.get_changed_df()
        self.param.set_parameters(param_df)
    
    def _set_callback(self):
        
        def run(change):
            nonlocal self
            print("running...")
            self.mySystem=CellSystem(parameters=self.param)
            with self._output:
                self.mySystem.run() 
        
        self._run_button.on_click(run)        
    
    def create_ui(self):
  
        return self.my_ui
    
#        self.solver.set_boundary("fixed")
#        self.solver.model="simple"
#        self.solver.N_r=600
#        t_min=1E-9
#        n_min=2000
#        t_max=3600*24*10
 #       n_step_pr_10=1000
 #       myT=Const_LogSpacing(t_min,t_max,n_step_pr_10,n_min)
 #       with self.output:
 #           print("Running....")
        #self.solver.solve(T=myT,profile={'type': 'soft','length':0.025})
        
    
    
#    def run(self):
#        self.solver=Glp_R_pde(parameters=self.param)
#        self.solver.set_boundary("fixed")
#        self.solver.model="simple"
#        self.solver.N_r=600
#        t_min=1E-9
#        n_min=2000
#        t_max=3600*24*10
 #       n_step_pr_10=1000
 #       myT=Const_LogSpacing(t_min,t_max,n_step_pr_10,n_min)
 #       with self.output:
 #           print("Running....")
        #self.solver.solve(T=myT,profile={'type': 'soft','length':0.025})
    

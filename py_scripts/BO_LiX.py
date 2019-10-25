from GPy.models import GPRegression
from GPy.kern.src.sde_matern import sde_Matern32 as Matern32
from GPy.kern.src.sde_matern import sde_Matern52 as Matern52
from GPy.kern.src.rbf import RBF
from emukit.model_wrappers import GPyModelWrapper
from emukit.core import ParameterSpace, ContinuousParameter
from emukit.experimental_design.model_free.latin_design import LatinDesign
from emukit.bayesian_optimization.loops import BayesianOptimizationLoop
from emukit.core.loop import UserFunctionWrapper
import numpy as np
from l2_bayes_opt.acquisitions import (
    L2NegativeLowerConfidenceBound as L2_LCB,
    L2ExpectedImprovement as L2_EI)
import matlab.engine
from Engine import LiX_wrapper

###

INIT_POINTS=20
BO_ITER = 20
NOISE = 0
frac_M = True
frac_X = True

###
eng = matlab.engine.start_matlab()

if NOISE == 0:
    NOISE = 1E-6 # numerical reason    

def save(bayes_state):
    X = bayes_state.loop_state.X
    Y = bayes_state.loop_state.Y
    
    return X,Y

def black_box_func(x): #!!! 
    global counter
    x = x.flatten()
    s_M = x[0]
    s_X = x[1]
    e_M = x[2]
    e_X = x[3]
    
    params = np.array([[s_M,s_X],[e_M,e_X]])
    
    #First True is the bo_mode data type
    opt_result = LiX_wrapper(True,'LiCl','Rocksalt','JC',params,
                             False,False,eng)
    
    counter += 1
    print('Iteration {}-->params:\n{}'.format(counter,
          np.array2string(params,precision=3)))
    print('Energy compared to target:{:.2f}'.format(opt_result[0][0]-target[0][0]))
    
    if np.isnan(opt_result[0][0]):
        opt_result = np.zeros((1,10))
    
    if frac_X == False:
        opt_result = opt_result[:,:7]
    if frac_M == False:
        opt_result = opt_result[:,:4]
    
    return opt_result
    
parameter_space = ParameterSpace([ContinuousParameter('sigma_M', 0.1, 0.5),
                                  ContinuousParameter('sigma_X', 0.1, 0.5),
                                  ContinuousParameter('epsilon_M', 1E-3, 1.5),
                                  ContinuousParameter('epsilon_X', 1E-3, 1.5)])
    
f=UserFunctionWrapper(black_box_func)

def main():
    print("######################")
    global target,X,Y0,values,frac_M,frac_X
    
    target_params = np.array([[0.14,0.4],[1.4,0.03]])
    
    # LiX_wrapper(bo_mode,salt,structure,model,JC_params,par,verbose,eng)
    target = LiX_wrapper(True,'LiCl','Rocksalt','JC',
                         target_params,False,False,eng)
    target = np.where(target==0,target+1,target)
    print('Target initialized!')
    
    latin_design = LatinDesign(parameter_space=parameter_space)
    X0 = latin_design.get_samples(INIT_POINTS)
    Y0 = np.array([])
    for x in X0:
        x = np.array([x])
        Y0 = np.append(Y0,f.evaluate(x))
    values = []
    for y in Y0:
        values.append(y.Y)
    
    values = np.asarray(values,dtype=float)
    values = values.reshape((-1,np.max(np.shape(target))))
    
    ### Redundancy check
    if (values[:,7:]==values[0,7]).all():
        values = values[:,:7]
        frac_X = False
        
    if (values[:,4:7]==values[0,4]).all():
        values = values[:,:4]
        frac_M = False
        
    ### BO Loop
    kern = Matern52(X0.shape[1],variance=1E-6)
    model = GPRegression(X0, values, kernel=kern,
                         normalizer=True, noise_var=NOISE) # Kernel = None: RBF default

    model.optimize(optimizer='lbfgsb')
    model.optimize_restarts(num_restarts=50,verbose=False)
    model_wrapped = GPyModelWrapper(model)
      
    acq = L2_LCB(model=model_wrapped, target=target, beta = np.float64(1)) 
    # beta is the exploration constant
    bayesopt_loop = BayesianOptimizationLoop(
                model=model_wrapped, space=parameter_space, acquisition=acq)
    bayesopt_loop.run_loop(f, BO_ITER)
    
    return save(bayesopt_loop)

if __name__ == '__main__':
    counter = 0
    test = main()
    eng.quit()
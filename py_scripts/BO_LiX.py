from GPy.models import GPRegression
from GPy.kern.src.sde_matern import sde_Matern32 as Matern32
from GPy.kern.src.sde_matern import sde_Matern52 as Matern52
from GPy.kern.src.rbf import RBF
from emukit.model_wrappers import GPyModelWrapper
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
salt = 'LiI'
structure = 'Rocksalt'
INIT_POINTS=30
BO_ITER = 30
NOISE = 1E-4
frac_M = True
frac_X = True
mix_rules = True
bo_flag = False
focus = None # None, 'energy ', or 'constant'
###
if not mix_rules:
    n_params = 6
else:
    n_params = 4
    
eng = matlab.engine.start_matlab()

if NOISE == 0:
    NOISE = 1E-6 # numerical reason    

def save(bayes_state):
    X = bayes_state.loop_state.X
    Y = bayes_state.loop_state.Y
    
    return X,Y

def black_box_func(x): #!!! 
    global counter, Y_original
    x = x.flatten()
    
    if max(np.shape(x)) == 4:
        s_M = x[0]
        s_X = x[1]
        e_M = x[2]
        e_X = x[3]
        
        params = np.array([[s_M,s_X],[e_M,e_X]])
        
    else:
        s_M = x[0]
        s_X = x[1]
        s_MX = x[2]
        e_M = x[3]
        e_X = x[4]
        e_MX = x[5]
        
        params = np.array([[s_M,s_X,s_MX],[e_M,e_X,e_MX]])
    
    #First True is the bo_mode data type
    # LiX_wrapper(bo_mode,salt,structure,model,model_params,par,verbose,eng)
    opt_result = LiX_wrapper(True,salt,structure,'JC',params,
                             False,False,eng)
    
    counter += 1
    Y_original = np.append(Y_original,opt_result)
    
    if np.isnan(opt_result[0][0]):
        opt_result = np.zeros((1,10))
        
    print('Iteration {}-->params:\n{}'.format(counter,
          np.array2string(params,precision=3)))
    print('Energy compared to target:{:.2f}'.format(opt_result[0][0]-target[0][0]))
    print('Lattice parameters compared to target:{:.2f}'.format(
                opt_result[0][1]-target[0][1]))
        
    
    if frac_X == False:
        opt_result = opt_result[:,:7]
    if frac_M == False:
        opt_result = opt_result[:,:4]
    
    if bo_flag:
        if focus == 'energy':
            return opt_result.flatten()[0].reshape(1,-1)
        if focus == 'constant':
            return opt_result.flatten()[1:4].reshape(1,-1)
        else:
            return opt_result[:,:].reshape(1,-1)
    else:
        return opt_result.reshape(1,-1)
    
if n_params == 4:
    from params_setting import parameter_space_4d as parameter_space
else:
    from params_setting import parameter_space_6d as parameter_space
    
f=UserFunctionWrapper(black_box_func)

def main():
    print("######################")
    global target,X0,Y0,values,frac_M,frac_X,bo_flag
    
    #target_params = np.array([[0.14,0.4],[1.4,0.03]])
    
    #target = LiX_wrapper(True,'LiF','Rocksalt','JC',
    #                     target_params,False,False,eng)

    target = np.array([[-764.5,6.012*0.99,6.012*0.99,6.012*0.99]])
    
    if focus == 'energy':
        target_comp = target[0,0].reshape(1,-1)
    if focus == 'constant':
        target_comp = target[0,1].reshape(1,-1)
    else:
        target_comp = target[0,:4].reshape(1,-1)
        
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

    ### Redundancy check
    if (values[:,7:-1]==values[0,7]).all():
        values = values[:,:7]
        frac_X = False
        
    if (values[:,4:7]==values[0,4]).all():
        values = values[:,:4]
        frac_M = False

    values = values.reshape(-1,np.max(np.shape(target)))    
    bo_flag = True
    
    if focus == 'energy':
        values = values[:,0].reshape(-1,1)
    if focus == 'constant':
        values = values[:,1:4].reshape(-1,3)
        
    ### BO Loop
    kern = Matern52(X0.shape[1],variance=1)
    model = GPRegression(X0, values, kernel=kern,
                         normalizer=True, noise_var=NOISE) # Kernel = None: RBF default

    model.optimize(optimizer='lbfgsb')
    model.optimize_restarts(num_restarts=50,verbose=False)
    model_wrapped = GPyModelWrapper(model)
      
    acq = L2_LCB(model=model_wrapped, target=target_comp, beta = np.float64(1.)) 
    # beta is the exploration constant
    bayesopt_loop = BayesianOptimizationLoop(
                model=model_wrapped, space=parameter_space, acquisition=acq)
    bayesopt_loop.run_loop(f, BO_ITER)
    
    return save(bayesopt_loop)

if __name__ == '__main__':
    counter = 0
    Y_original = np.array([])
    test = main()
    Y_original = Y_original.reshape(INIT_POINTS+BO_ITER,-1)
    
    tmp = (Y_original,test[0],test[1],'LiI','Rocksalt','with_mix','both','30+30')
    import pickle
    with open('./tmp_data/LiI/opt_both_MIX.pickle', 'wb') as f:
        pickle.dump(tmp, f)
    eng.quit()
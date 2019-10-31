import numpy as np
import matlab.engine
import io

def np2mat(np_array):
    if np_array.ndim != 2:
        raise ValueError('Expected 2d array, but got {}d array.'.format(
                np_array.ndim))
    tmp_np = np_array.tolist()
    
    return matlab.double(tmp_np)

def LiX_mat2py(salt,structure,model,JC_params,par=False,verbose=False,eng=None):
    eng.cd(r'..')
    
    out = io.StringIO()
    err = io.StringIO()
    
    if eng == None:
        raise RuntimeError('Please pass the matlab engine.')
        
    if not verbose:    
        if not par:
            data = eng.Structure_Minimization(salt,structure,model,
                                              np2mat(JC_params),True,
                                              nargout=1,
                                              stdout=out,stderr=err)
        else:
            data = eng.Structure_Minimization_Par(salt,model,
                                                  np2mat(JC_params),True,
                                                  nargout=1,
                                                  stdout=out,stderr=err)
    else:
        if not par:
            data = eng.Structure_Minimization(salt,structure,model,
                                              np2mat(JC_params),True,
                                              nargout=1)
        else:
            data = eng.Structure_Minimization_Par(salt,model,
                                                  np2mat(JC_params),True,
                                                  nargout=1)
    eng.cd(r'./py_scripts')
    
    return np.asarray(data,dtype='float64')

def LiX_wrapper(bo_mode,*args):
    tmp_data = LiX_mat2py(*args)
    
    if not bo_mode:
        result = {
                'energy':tmp_data[0][0],
                'latt_params':tmp_data[0][1:4],
                'frac_coord_M':tmp_data[0][4:7],
                'frac_coord_X':tmp_data[0][7:]
                }
        return result
    else:
        return tmp_data
    
eng = matlab.engine.start_matlab()
params = np.array([[0.11,0.17],[0.038,1.088]])
test = LiX_wrapper(True,'LiCl','Rocksalt','JC',params,False,True,eng)
eng.quit()

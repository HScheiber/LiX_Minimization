import numpy as np
import matlab.engine
import io

def np2mat(np_array):
    if np.shape(np_array) not in [(2,2),(2,3)]:
        np_array = np_array.reshape(2,-1)
    tmp_np = np_array.tolist()
    
    return matlab.double(tmp_np)

def trans_2d(a,n_points,n_params,debug=False):
    if n_params <10:
        a = a.reshape(n_points,-1)
        n_params = int(n_params/2)
        tmp = []
        for array in a:
            tmp.append(array.reshape(-1,n_params).T)
        
        tmp = np.asarray(tmp).T
        if not debug:
            return matlab.double(tmp.tolist())
        else:
            return tmp
    else:
        a = a.reshape(n_points,-1)
        n_params = 3
        tmp = []
        for array in a:
            tmp.append(array.reshape(-1,n_params).T)
        
        tmp = np.asarray(tmp).T
        if not debug:
            return matlab.double(tmp.tolist())
        else:
            return tmp
        
def LiX_mat2py(salt,structure,model,params,par=False,verbose=False,eng=None,
               inner_opt=False,CRDamping=True,C6Damping=0,background=False):
    eng.cd(r'..')
    
    out = io.StringIO()
    err = io.StringIO()
    
    if eng == None:
        raise RuntimeError('Please pass the matlab engine.')
        
    if not verbose:    
        if not par:
            data = eng.Structure_Minimization(salt,structure,model,
                                              np2mat(params),inner_opt,CRDamping,
                                              C6Damping,nargout=1,background=background,
                                              stdout=out,stderr=err)
        else:
            data = eng.Structure_Minimization_Par(salt,model,
                                                  np2mat(params),inner_opt,CRDamping,
                                                  C6Damping,nargout=1,background=background,
                                                  stdout=out,stderr=err)
    else:
        if not par:
            data = eng.Structure_Minimization(salt,structure,model,
                                              np2mat(params),inner_opt,CRDamping,
                                              C6Damping,nargout=1,background=background)
        else:
            data = eng.Structure_Minimization_Par(salt,model,
                                                  np2mat(params),inner_opt,CRDamping,
                                                  C6Damping,nargout=1,background=background)
    eng.cd(r'./py_scripts')
    
    return np.asarray(data,dtype='float64')

def LiX_mat2py_multi(salt,structure,model,parameters,CRDamping=True,
                     C6Damping=0,verbose=False,eng=None):
    eng.cd(r'..')
    
    out = io.StringIO()
    err = io.StringIO()
    
    if eng == None:
        raise RuntimeError('Please pass the matlab engine.')
    
    n_points,n_params = parameters.shape
    
    if not verbose:
        data = eng.Structure_Minimization_Multi(salt,structure,model,
                                               trans_2d(parameters,n_points,n_params),False,CRDamping,
                                               C6Damping,nargout=1,
                                               stdout=out,stderr=err)
    else:
        data = eng.Structure_Minimization_Multi(salt,structure,model,
                                              trans_2d(parameters,n_points,n_params),False,CRDamping,
                                              C6Damping,nargout=1)
    eng.cd(r'./py_scripts')
    
    return np.asarray(data,dtype='float64')
    
def LiX_wrapper(bo_mode,*args,**kwargs):
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

'''    
eng = matlab.engine.start_matlab()
params = np.array([[0.2250,0.4508,0.2966,1.16,0.4706,0.7086]]).reshape(2,3)
test = LiX_wrapper(True,'LiI','Rocksalt','JC',params,False,True,eng)
eng.quit()
'''
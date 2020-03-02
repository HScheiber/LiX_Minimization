import numpy as np
from Engine import LiX_wrapper
import matlab.engine
from Engine import np2mat

#LiX_wrapper(bo_mode,salt,structure,model,JC_params,par,verbose,eng)
if '__main__' == __name__:
    eng = matlab.engine.start_matlab()
    X0 = np.load('./X0.npy')
    eng.cd('..')
    par_test = eng.Structure_Minimization_Multi('LiF','Rocksalt','JC',
                                                np2mat(X0),True,nargout=1)
    eng.quit()
    
    np.save('./par_result.npy',par_test)
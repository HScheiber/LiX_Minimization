import numpy as np
from Engine import LiX_wrapper
import matlab.engine

#LiX_wrapper(bo_mode,salt,structure,model,JC_params,par,verbose,eng)
eng = matlab.engine.start_matlab()
params = np.array([[0.14,0.4],[1.4,0.03]])
par_test = LiX_wrapper(True,'LiCl','Rocksalt','JC',params,True,True,eng)
eng.quit()

np.save('./par_result.npy',par_test)
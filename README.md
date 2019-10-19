# LiX_Minimization
 Steepest descent minimization script which interfaces with GROMACS to produce optimized LiX crystal structures.

All scripts and files necessary to run the primary minimization script are included in the package.
The main script is called “Structure_Minimization.m”. This takes as input: (Salts,Structures,Models,Scale_Dispersion,Scale_Epsilon,Glost_On,OptPos)
Salts is a cell array of strings, indicating which lithium halide salts to optimize. Options are ‘LiF’ ‘LiCl’ ‘LiBr’ ‘LiI’.
Structures is a cell array of strings, indicating which crystal structures to optimize. Options are 'Rocksalt' 'Wurtzite' 'Sphalerite' 'CsCl' 'NiAs' 'BetaBeO' 'FiveFive'.
Models is a cell array of strings, indicating which models to use. Current options are ‘JC’ (which is JC model for SPC/E water), ‘JC3P’ (which is the JC model for TIP3P water), ‘JC4P’ (which is the JC model for TIP4P_EW water), and ‘TF’ (which is the Tosi-Fumi model).
	In the future, I think it will be most reasonable to modify this part of the input to be an array of model parameters.
Scale_Dispersion is an array of floats, indicating how much to scale the dispersion parameters of a particular model. This should be removed for your project.
Scale_Epsilon is an array of floats, indicating how much to scale the epsilon parameter of a model (only works for JC models). This should be removed for your project.
Glost_On is a Boolean. If TRUE, it forces the script to run GROMACS on a single thread. I used this for running parallel jobs with GLOST (Greedy Launcher of Small Tasks: https://docs.computecanada.ca/wiki/GLOST). 
GLOST is less buggy than Matlab’s inbuilt parallel mode when running on compute Canada clusters, because certain aspects of Matlab parallelization cannot be modified by regular users, making it very annoying to run properly.
When this option is set to FALSE, the GROMACS software will try to run on as many cores as are available.
OptPos is a Boolean. When TRUE it allows for the optimization of lattice parameters AND fractional coordinates. When FALSE it only optimizes lattice parameters, fixing the fractional coordinates in place.
	Note that when fractional coordinates are fixed, one structure will never converge towards a different structure. Hence this option was useful if I wanted to get (for example) an idea of LiF wurtzite energy, even though it wasn’t a true local minimum. 
	Also, this script maintains the overall symmetry of the space group. Therefore, all unit cell angles are held fixed for all structures, and some relationships between cell parameters are fixed as well (e.g. rocksalt a = b = c).
	It would not be too much work for me to remove this constraint, but I very much doubt it would change the results significantly. 
I’ve done some ab initio calculations of these structures in space group P1 (i.e. all cell angles and cell parameters are free to change) and structures did not change significantly.

Currently the output of this script is just what’s printed out when you run it. 

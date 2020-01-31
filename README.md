# LiX_Minimization
Steepest descent minimization script which interfaces with GROMACS to produce optimized LiX crystal structures.

All scripts and files necessary to run the primary minimization script are included in the package.

The main script is called “Structure_Minimization.m�?. This takes as input:
`(Salt,Structure,Model,Parameters,OptPos)`


# INFO ABOUT INPUTS
`Model` is a string. One of 'JC' or 'TF'. Sets the mathematical form of the model to use.

`Salt` is a string. Can be any of: 'LiF' 'LiCl' 'LiBr' 'LiI' 'NaCl'

`Structure` is a string. Can be any of: 'Rocksalt' 'Wurtzite' 'Sphalerite' 'CsCl' 'NiAs' 'BetaBeO' 'FiveFive'

`OptPos` is a Boolean. When TRUE it allows for the optimization of lattice parameters AND fractional coordinates. When FALSE it only optimizes lattice parameters, fixing the fractional coordinates in place.
* Note that when fractional coordinates are fixed, one structure will never converge towards a different structure. Hence this option was useful if I wanted to get (for example) an idea of LiF wurtzite energy, even though it wasn’t a true local minimum. 
* Also, this script maintains the overall symmetry of the space group. Therefore, all unit cell angles are held fixed for all structures, and some relationships between cell parameters are fixed as well (e.g. rocksalt a = b = c).
* It would not be too much work for me to remove this constraint, but I very much doubt it would change the results significantly. 


`Parameters` is and `N x M` matrix of floats, `N` and `M` depend on the chosen model. 

For JC model: Parameters are contained in either a `2 x 2` OR `2 x 3` array.
For TF model: Parameters are contained in a `4 x 3` array.

Below are descriptions of the Parameter matrix options. M refers to the
Metal ion (i.e. Li or Na) and X refers to the Halide ion (i.e. F, Cl, Br,
or I)

JC 2 x 2 matrix:  
<pre>
σ_M     σ_X        (units: nm)  
ε_M     ε_X        (units: kJ/mol)
</pre>

When using a 2 x 2 array with the JC model, 
the Lorenz-Berthelot mixing rules are assumed.


JC 2 x 3 array:  
<pre>
σ_MM    σ_XX    σ_MX       (units: nm)  
ε_MM    ε_XX    ε_MX       (units: kJ/mol)
</pre>

When using a 2 x 3 array with the JC model, mixing rules are not used.


TF 4 x 3 array:  
<pre>
α_MM    α_XX    α_MX       (units: nm^-1)  
B_MM    B_XX    B_MX       (units: kJ/mol)  
C_MM    C_XX    C_MX       (units: (kJ nm^6)/mol)  
D_MM    D_XX    D_MX       (units: (kJ nm^8)/mol)
</pre>

No mixing rules are defined for the TF model.

`CRDamping` Is a boolean that, when true, switches on a close-range logistic damping function on all attractive interactions (including opposite sign coulombic interactions). This prevents the code from creating unphysically tightly bound systems, which can cause the program run slow or crash. The damping function is designed not to affect the pair potential at physically relavent distances.

The output from the Structure_Minimization script is a 1 x 10 numeric array which contains the following information in order:  
`Output[1]  = Lattice energy in kJ/mol of formula units (i.e. energy per mole of ion pairs)`  
`Output[2]  = Length of parameter a`  
`Output[3]  = Length of parameter b`  
`Output[4]  = Length of parameter c`  
`Output[5]  = Fractional coordinate x of the metal ion within the asymmetric unit`  
`Output[6]  = Fractional coordinate y of the metal ion within the asymmetric unit`  
`Output[7]  = Fractional coordinate z of the metal ion within the asymmetric unit`  
`Output[8]  = Fractional coordinate x of the halide ion within the asymmetric unit`  
`Output[9]  = Fractional coordinate y of the halide ion within the asymmetric unit`  
`Output[10] = Fractional coordinate z of the halide ion within the asymmetric unit`  

This information (plus the input structure) is enough to uniquely determine the unit cell of any of the 7 candidate structures, as long as the space group of the unit cell is held fixed.
It turns out that all 7 of the candidate structures have only a single ion pair within their asymmetric unit.

I’ve done some ab initio calculations of these LiX structures in space group P1 (i.e. all cell angles and cell parameters are free to change) and optimized structures did not change significantly.

Currently the output of this script is just what’s printed out when you run it. There is also data files that are produced, but I have set the script to delete these at the end. When I was collecting data I used to have a separate script that would scrape the data files and save them as matlab data structures.

In the subdirectory called `DATA` is all of this previous data I've collected. This will be useful in generating reasonable initial conditions.

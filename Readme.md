# Boussinesq thermal dynamo in a full sphere

The model equations are

![Model equations](https://latex.codecogs.com/gif.latex?%5Cinline%20%5Cbegin%7Balign*%7D%20%5Cleft%28%5Cpartial_t%20-%20Pm%20%5Cboldsymbol%7B%5CDelta%7D%5Cright%29%5Cmathbf%7Bu%7D%20%26%20%3D%20%5Cmathbf%7Bu%7D%20%5Ctimes%20%5Cleft%28%20%5Cboldsymbol%7B%5Cnabla%7D%20%5Ctimes%20%5Cmathbf%7Bu%7D%5Cright%29%20&plus;%20%5Cfrac%7BPm%5E2%20Ra%7D%7BE%20Pr%7D%20%5CTheta%20%5Cmathbf%7Br%7D%20-%5Cfrac%7BPm%7D%7BE%7D%5Chat%7B%5Cmathbf%7Bz%7D%7D%5Ctimes%5Cmathbf%7Bu%7D%20&plus;%20%5Cfrac%7BPm%7D%7BE%7D%5Cleft%28%5Cboldsymbol%7B%5Cnabla%7D%5Ctimes%5Cmathbf%7BB%7D%5Cright%29%5Ctimes%20%5Cmathbf%7BB%7D%20-%20%5Cnabla%5CPi%20%5C%5C%5B0.6cm%5D%20%5Cleft%28%5Cpartial_t%20-%20%5Cboldsymbol%7B%5CDelta%7D%5Cright%29%5Cmathbf%7BB%7D%20%26%20%3D%20%5Cboldsymbol%7B%5Cnabla%7D%20%5Ctimes%20%5Cleft%28%5Cmathbf%7Bu%7D%20%5Ctimes%20%5Cmathbf%7BB%7D%5Cright%29%5C%5C%5B0.6cm%5D%20%5Cleft%28%5Cpartial_t%20-%20%5Cfrac%7BPm%7D%7BPr%7D%5CDelta%5Cright%29%5CTheta%20%26%20%3D%20S%20-%20%5Cmathbf%7Bu%7D%5Ccdot%5Cnabla%5CTheta%5C%5C%5B0.6cm%5D%20%5Cboldsymbol%7B%5Cnabla%7D%5Ccdot%5Cmathbf%7Bu%7D%20%26%20%3D%200%5C%5C%5B0.6cm%5D%20%5Cboldsymbol%7B%5Cnabla%7D%5Ccdot%5Cmathbf%7BB%7D%20%26%20%3D%200%20%5Cend%7Balign*%7D)

with the parameters defined as

![Nondimensional parameters](https://latex.codecogs.com/gif.latex?%5Cinline%20%5Cbegin%7Balign*%7D%20Pr%20%26%20%3D%20%5Cfrac%7B%5Cnu%7D%7B%5Ckappa%7D%5C%5C%20Pm%20%26%20%3D%20%5Cfrac%7B%5Cnu%7D%7B%5Ceta%7D%5C%5C%20Ra%20%26%20%3D%20%5Cfrac%7Bg%20%5Calpha%20%5Cbeta%20r_o%5E4%7D%7B2%5COmega%5Ckappa%7D%5C%5C%20E%20%26%20%3D%20%5Cfrac%7B%5Cnu%7D%7B2%5COmega%20r_0%5E2%7D%20%5Cend%7Balign*%7D)

The system is solved using a Toroidal/Poloidal decomposition of the velocity ![u](https://latex.codecogs.com/svg.latex?%5Cinline%20%5Cmathbf%7Bu%7D):

![Toroidal/Poloidal decomposition](https://latex.codecogs.com/svg.latex?%5Cinline%20%5Cmathbf%7Bu%7D%3D%5Cmathbf%7B%5Cnabla%7D%5Ctimes%20T%20%5Cmathbf%7Br%7D%20&plus;%20%5Cmathbf%7B%5Cnabla%7D%5Ctimes%5Cmathbf%7B%5Cnabla%7D%5Ctimes%20P%20%5Cmathbf%7Br%7D)

and for the magnetic field ![B](https://latex.codecogs.com/svg.latex?%5Cinline%20%5Cmathbf%7BB%7D):

![Toroidal/Poloidal decomposition](https://latex.codecogs.com/gif.latex?%5Cinline%20%5Cmathbf%7BB%7D%3D%5Cmathbf%7B%5Cnabla%7D%5Ctimes%20%5Cmathcal%7BT%7D%20%5Cmathbf%7Br%7D%20&plus;%20%5Cmathbf%7B%5Cnabla%7D%5Ctimes%5Cmathbf%7B%5Cnabla%7D%5Ctimes%20%5Cmathcal%7BP%7D%20%5Cmathbf%7Br%7D)

Note: To access the equations in codecogs online LaTeX editor replace "https://latex.codecogs.com/svg.latex?" with "https://www.codecogs.com/eqnedit.php?latex=".

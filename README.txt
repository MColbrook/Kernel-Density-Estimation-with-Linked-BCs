This folder contains code for the paper:

KERNEL DENSITY ESTIMATION WITH LINKED BOUNDARY CONDITIONS

authored by Matthew Colbrook, Zdravko Botev, Karsten Kuritz and Shev MacNamara. 

See license.txt.

Description for the files is as follows.

Cts_samples: 
For most users, THIS IS THE FUNCTION THAT YOU WANT. 
This function will estimate a density, given your data.
This is similar to Cts_func but here the input is a vector D of samples.



Cts_func: 
given input (f0,T,r,x,PLOT), this solves for the continuous model at the points x.
f0 is the initial values (function handle), T the stopping time, r the BCs
parameter and PLOT the option of plotting the result



solution_k: 
summation function used by the continuous model. 
This is an internal function (that helps some of the other functions and so
is required to be in the same directory). Most users will never want to
directly open or call this function.

Discrete_linked: 
Similar to Cts_samples but solves using the discrete matrix approach.
There is an additional parameter m which determines the number of (equally
spaced) bins for the data.
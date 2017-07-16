# L-shaped-LDC-
Code to compute flows in a single-sided L-shaped cavity Using Lattice Boltzmann Technique (TRT) model
using a D2Q9 model

Boundary conditions used:
Zou-He boundary condition on top moving wall
and Bounce-back on remaining stationary walls

Geometry of L-shaped cavity can be found in L-cavity_geo.pdf file in this repository.
The aspect ratio is defined as : ![equation](http://www.sciweavers.org/tex2img.php?eq=AR%20%3D%20%5Cfrac%7BL%20-%20L_1%7D%7BL%7D&bc=White&fc=Black&im=jpg&fs=12&ff=arev&edit=0), according to the nomenclature as shown in the geometry. 

Code modified from Jonas Latt's code cavity2d.m
available at http://wiki.palabos.org/numerics:matlab_samples

More details of TRT-LBM can be found in the [Article by Irina Ginzburg](https://www.researchgate.net/profile/Irina_Ginzburg/publication/281975432_Two-relaxation-time_Lattice_Boltzmann_scheme_about_parametrization_velocity_pressure_and_mixed_boundary_conditions/links/5600745d08ae07629e52adc0/Two-relaxation-time-Lattice-Boltzmann-scheme-about-parametrization-velocity-pressure-and-mixed-boundary-conditions.pdf)

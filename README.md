# Two_dielectric_sphere_electrostatic_force_by_normal_deri
A matlab script for post-processing the BEM output; give the force on two dielectric spheres with central charges.
Two input files,
data.dat, contains the BEM outputs on the collocation points, it includes the location of the pts, the normal derivative on these pts, and the weights for each pts (if we use CC, then the weights are the triangle areas).
source.dat, contains the information about the two spheres, including the location of their centers, their central charges, and the inside/outside dielectric constants.

The output gives the force on both spheres, in x, y, z directions. It is compared with Hybrid method; and for the current test files (calculated using CC and 1280 triangles on each sphere), the accuracy in force is more than 4 digits.

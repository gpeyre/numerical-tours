Toolbox Wavelets on Meshes - wavelet transform on 3D meshes

This toolbox allows to compute the wavelet transform of a function defined on semi-regular triangulation. For instance it allows to compute
* wavelet transform of a function defined on the sphere, one this function has been sampled on a 4:1 subdivided tetrahedron.
* wavelet transform of a function defined on a 4:1 subdivided coarse triangulation.
* wavelet transform of a semi-regula meshes, viewed as 3 functions defined on a 4:1 subdivided coarse triangulation.

The wavelet transform is implemented using the lifting scheme, as described in 
       	Peter Schrodder and Wim Sweldens
 		Spherical Wavelets: Efficiently Representing Functions on the Sphere
       	Siggraph 95

This toolbox also allows to compute a semi-regular triangular subdivision surface using Loop and Butterfly stencils.

There are helper functions to create multiresolution spherical meshes, to load a semi-regular mesh from a geometry image, and to display function on sphere and on meshes.

This toolbox is still experimental. Only a lifted wavelet using a butterfly sub-division is implemented - more to come soon (including Loop wavelets). 

Copyright (c) 2007 Gabriel Peyre
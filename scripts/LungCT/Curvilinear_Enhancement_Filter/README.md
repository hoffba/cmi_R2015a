
README.txt  (September 03, 2012)
--------------------------

The contents of the rock-segmentation.zip file gets extracted into a folder called "RockSeg".

# Curvilinear Enhancement Filter
3D and 2D Hessian based tubular (vessel/vesselness) and spherical (blob/blobness) enhancement filters.

This code is a Matlab implementation of the algorithm described in the Journal paper by Sundaresh Ram and Jeffrey J. Rodriguez below.

1. [Ram S, Danford F, Howerton S, Rodriguez JJ, Geest JPV., "*Three-Dimensional Segmentation of the Ex-Vivo Anterior Lamina Cribrosa From Second-Harmonic Imaging Microscopy*", IEEE Transactions on Biomedical Engineering, 65(7), p. 1617-1629 (2018), doi={10.1109/TBME.2017.2674521}]


The code is based on Dirk-Jan Kroon's implementation of Frangi's vesselness filter. (https://www.mathworks.com/matlabcentral/fileexchange/24409-hessian-based-frangi-vesselness-filter)

### Tips:

* Make sure that the objects of interest have the highest (if bright compared to the background) or lowest (if dark compared to background) intensities in the image/volume. Scale/normalize the images appropriately.

* The 3D method contains a c-code file that needs to be compiled with "mex eig3volume.c". (For more info visit: https://www.mathworks.com/matlabcentral/fileexchange/24409-hessian-based-frangi-vesselness-filter)

* Threshold the filter response to remove any remaining enhanced noise

### Content:

2D enhancement of vessel/tube-like structures:

 * vesselness2D.m - main function
 * example_vesselness2D.m - filter applied on a 2D retinal vasculature
 * fundus2D.png - image for the example
 
3D enhancement of vessel/tube-like structures:

 * vesselness3D.m - main function
 * eig3volume.c - fast computation of eigenvalues
 * example_vesselness3D.m - filter applied on a 3D cerebral vasculature
 * volume.mat - volume for the example
 
3D enhancement of blob/sphere-like structures:

 * blobness3D.m - main function
 * eig3volume.c - (as above)
 * example_blobness3D.m - aneurysm enhancement in a 3D cerebral vasculature
 * volume.mat - (as above)
 
 2D enhancement of blob/disk-like structures (**experimental - not published**):

 * blobness2D.m - main function
 * example_blobness2D.m - enhancement of an aneurysm on a vessel (*synthetically generated image*)
 * blob2D.png - image for the example
  

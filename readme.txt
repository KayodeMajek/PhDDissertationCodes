Written by Olukayode T. Majekodunmi for PhD Dissertation titled "Discontinuous Colloidal Clogging in Tapered Microchannels" submitted in May, 2024. 

"ImageAnalyzer.m" works alongside the function m-files "sppcrop.m" and "nonclogs.m"; they should be downloaded, stored and run from the same folder. The code measures the dimensions of clogs formed after a monodisperse suspension of fluorescent particles is flowed through a system of parallel microchannels. The channels are tapered in this application. 

1. The code collects all image files of the type "TIFF" from a directory: f = 'C:\Users\Tilescan\*.tif'. The directory should be changed as needed.
2. The images are first tiled. The code identifies all the clogs, which are 8-connected objects in the tiled image. The location and dimensions of the clogs are measured and stored.

3. "sppcrop.m" reruns the "regionprops" function to identify the clogs. This is necessary as empty/blank elements in the end positions of tiled image array are deleted in the analysis. These are blank spaces beyond the outlet of the channels, empty top and bottom sections (no channels included) of the microfluidic device. This adjustment does not change the clog dimensions. However, they change the locations in the final analysis. The purpose here is to have a form of uniformity across the different conditions studied.

4. A clog should be the same size as the width of the channel corresponding to its location. However, there are instances where particles not associated with any clog are identified as clogs. "nonclogs.m" removes information about these objects from the collected data. These objects are often particles not associated with any clog in the channels as their width as much smaller than the channel width corresponding to their location.

For more information, comments, and potential collaborations, please contact:
Prof. Sara Hashmi,
Complex Fluids Lab, Department of Chemical Engineering, Northeastern University. 
s.hashmi@northeastern.edu || https://github.com/smhashmi786/hashmilab || https://hashmilab.sites.northeastern.edu/

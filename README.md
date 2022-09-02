# EBE_Extended_Biomass_Estimation
Matlab tool for computing the biomass of riparian vegetation based on georeferenced point clouds (e.g., LiDAR).
_______________________________________________________________________________________________________________
 Extended Biomass Estimation (EBE) - version 2.0
  
  Matlab (created with Matlab 2021b)
 
  Melissa Latella & Tommaso Raimondo
 
  Created: November 2021
  Last update: 26 July 2022
 
  melissa.latella@polito.it
  
  Cite as: "Latella, M., Raimondo, T., Belcore, L., Salerno, L., 
                        and Camporeale, C. (2022), On the integration of 
                        LiDAR and field data for riparian biomass 
                        estimation, Journal of Environmental Management".
  
  and "Latella, M., Sola, F., & Camporeale, C. (2021). A Density-Based 
           Algorithm for the Detection of Individual Trees from LiDAR Data. 
            Remote Sensing, 13(2), 322".
 
  Aim: Mapping riparian vegetation biomass considering trees and shrubs along river
        channel(s) and floodplain by combining field data and LiDAR point clouds.
  
  Modules: 
  1. EBE.m is the module to map trees, shurbs and bushes along river channels 
  and floodplains. Based on "Latella, M., Raimondo, T., Belcore, L., Salerno, L., 
  and Camporeale, C. (2022), On the integration of LiDAR and field data for riparian 
  biomass estimation".
  2. EBE_ITDM.m is the module for individual tree detection and measurement
 	based on "Latella, M., Sola, F., & Camporeale, C. (2021). A 
 	Density-Based Algorithm for the Detection of Individual Trees from 
  LiDAR Data. Remote Sensing, 13(2), 322.

  Auxiliary functions:
  3. EBE_concave_hull.m is an auxiliary function for the EBE.m module. It is based 
  on "Moreira, A., & Santos, M. Y. (2007). Concave hull: A k-nearest neighbours 
 	approach for the computation of the region occupied by a set of points, 
  GRAPP 2007 - International Conference on Computer Graphics Theory and Application"
  and "Yan, Z., Liu, R., Cheng, L., Zhou, X., Ruan, X., & Xiao, Y. (2019). 
  A concave hull methodology for calculating the crown volume of individual trees 
  based on vehicle-borne LiDAR data. Remote Sensing, 11(6), 623."
  4. jsonInitialize.m is a function to initialize the json output file of the
  EBE.m module
  5. jsonUpdate.m is a function to update the json output file of the  EBE.m module
  
  Dependencies (EBE.m module):
  - Parallel Toolbox™ by Matlab (not mandatory)
  - Lidar Toolbox™ by Matlab (mandatory)
  - Individual tree position known or output from EBE_ITDM.m (not mandatory).
  - jsonInitialize.m and jsonUpdate.m (mandatory)
  - EBE_concave_hull.m (mandatory). In turn, the function demands selfintersect.m 
    by Antoni J. Canós (2022). Fast and Robust Self-Intersections 
    (https://www.mathworks.com/matlabcentral/fileexchange/13351-fast-and-robust-self-intersections), 
    MATLAB Central File Exchange. Retrieved July 26, 2022.
 
  Dependencies (ITDM.m module):
  - Parallel Toolbox™ by Matlab (not mandatory)
 
  The "Example" provides the input and output of a toy version of the model. Please 
  notice that the input files "Input/RGB/1.tif" and "Input/pointClouds/1.las" were
  not uploaded on GitHub due to their size, but can be found on Zenodo.
  

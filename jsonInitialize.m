%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Extended Biomass Estimation EBE - version 2.0
%
% Module to initialize the json output file
%
% Melissa Latella
%
% Created: 01 November 2021
% Last update: 26 July 2022
%
% melissa.latella@polito.it
%
% Cite as: "Latella, M., Raimondo, T., Belcore, L., Salerno, L.,
%                       and Camporeale, C. (2022), On the integration of
%                       LiDAR and field data for riparian biomass
%                       estimation, Journal of Environmental Management".
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function jsonData = jsonInitialize(jsonName,crs)

    clear jsonData
    jsonData.type='FeatureCollection';
    jsonData.name= jsonName;
    jsonData.crs=[];
    jsonData.crs.type="name";
    jsonData.crs.properties=[];
    jsonData.crs.properties.name= strcat("urn:ogc:def:crs:EPSG::",string(crs));
    jsonData.features=[];
    jsonData.features ={};

end
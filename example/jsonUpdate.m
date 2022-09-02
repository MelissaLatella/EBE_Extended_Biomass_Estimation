%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Extended Biomass Estimation EBE - version 2.0
%
% Module to update the json output file
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

function data = jsonUpdate(data, vegType, id, biomass, polygon)
           
            %single plant
            data0.features=[];
            data0.features.type = "Feature";
            data0.features.properties.veg =  vegType;
            data0.features.properties.id =  string(id);
            data0.features.properties.biomass = biomass;
            data0.features.geometry.type = "MultiPolygon";
            data0.features.geometry.coordinates = [polygon.Vertices(:,1) polygon.Vertices(:,2)];
            data0.features.geometry.coordinates = {data0.features.geometry.coordinates};
            data0.features.geometry.coordinates = {data0.features.geometry.coordinates};
            
            %array of plants
            data.features = [data.features;data0.features];
end

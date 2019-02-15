function [areas,colorMap,abbreviations]=visualAreas(addV1)
% Returns the visual areas of interest as a cell array that corresponds to names in the Allen Atlas 
%
% function [areas,colorMap,abbreviations] = brainAreaNames.visualAreas
%
%
% Purpose
% This is just a simple way of returning the list of area names we care about.
%
%
% Inputs
% addV1 - true by default. Shades V1 in gray.
% 
% Usage example:
% [areas,colorMap,abbreviations] = brainAreaNames.visualAreas(addV1)
% 
%
% Rob Campbell - Basel 2017

if nargin==0
    addV1=true;
end


areas = {...
    'Anteromedial visual area',...
    'Anterolateral visual area',...
    'posteromedial visual area',...
    'Lateral visual area',... %LM
    'Temporal association areas',...
    'Anterior area',...
    'Posterolateral visual area',...
    'Postrhinal area',... 
    'Perirhinal area',...
    'Rostrolateral visual area',...
    'Ectorhinal area',...
    'Laterointermediate area',...
    };

abbreviations = {...
    'AM',...
    'AL',...
    'PM',...
    'LM ',.... 
    'TEA',...
    'A ',...
    'P',... 
    'POR',... 
    'PER',... %OK?
    'RL',... 
    'ECT',...
    'LI',...
    'RS',...
    };


if addV1
    abbreviations=['V1',abbreviations];
    areas=['Primary visual area',areas];
end


%Create color tables for these areas and cells which project to them
n=length(areas);

%Area colours
myMap = [1,0.85,0.85;... %AM is pastel red
        0.85,1.0,0.85;... %AL is pastel green
        0.7,1.0,1.0;...  %PM is cyan
        1.0,0.66,1.0;...  %LM is magenta
        1.0,0.75,0.25;... %TEA is orange
        1.0,0.95,0.63;...  %Make the anterior area yellow
        0.0,1.0,0.0;... %posteriolateral is bright green
        1,0.5,0.5]; %post-rhinal is brighter red


if addV1
    myMap=[0.92,0.92,0.92; myMap]; %V1 is gray (addV1 must be true)
end

myMap = [myMap; colormap(brighten(parula(n-size(myMap,1)),0.85)) ];


colorMapAreas= myMap;
colorMapCells = colormap(brighten(myMap,-0.5));



colorMap.areas = containers.Map;
colorMap.cells = containers.Map;

%To reference using ID rather than name
SL=getAllenStructureList;

colorMap.areasID = containers.Map('KeyType','double','ValueType','any');
colorMap.cellsID = containers.Map('KeyType','double','ValueType','any');

for ii=1:n
    areaName = strrep(areas{ii},' ','_');
    colorMap.areas(areaName) = colorMapAreas(ii,:);
    colorMap.cells(areaName) = colorMapCells(ii,:);

    f=strmatch(areas{ii},SL.name,'exact');
    if isempty(f)
        fprintf('%s failed to find area %s. Skipping. \n', areas{ii});
        continue
    end

    colorMap.areasID(SL.id(f)) = colorMapAreas(ii,:);
    colorMap.cellsID(SL.id(f)) = colorMapCells(ii,:);
end




colorMap.areaNames=areas;

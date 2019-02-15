function [pltData,areasProcessed] = filterPlotAreasFromPltData(pltData,areas,retain)
% remove or select the listed brain areas from the pltData structure 
%
% function [pltData,excluded] = filterPlotAreasFromPltData(pltData,areas,retain)
%
% Returns a cell array of strings with the exluded areas and the 
% data with the removed areas.
%
% if retain is true (false by default) then we keep only "areas". Otherwise we
% remove "areas"
%
%
%
% Rob Campbell
%
% see - pointsByAreaPlot, clusterTracedCells

if nargin<3
    retain=0;
end

areasProcessed = {}; %list of areas we managed to exclude
if isempty(areas)
    return 
end

if ~iscell(areas)
    fprintf('Argument excludeAreas should be a cell array. Ignoring it.\n')
    return
end


if retain
    areasToKeep=[];
end

areasProcessed={};
verbose=1;

if verbose
    fprintf('\n\nEntering %s\n', mfilename)
end

for ii=1:length(areas)
    thisArea = areas{ii};

    ind=strmatch(thisArea,pltData.areaNamesInSamples);

    if isempty(ind)
        fprintf('Could not find area matching string %s. Skipping\n',thisArea)
        continue
    else
        if verbose
            fprintf('Keeping %s\n',pltData.areaNamesInSamples{ind})
        end     
    end

    areasProcessed = [areasProcessed; pltData.areaNamesInSamples(ind)];

    if ~retain %we exclude data
        pltData.allAreas(ind,:)=[];
        pltData.dataMat(ind,:)=[];
        pltData.areaNamesInSamples(ind)=[];
    else
        areasToKeep = [areasToKeep; ind];
    end

end

if retain
    pltData.allAreas=pltData.allAreas(areasToKeep,:);
    pltData.dataMat=pltData.dataMat(areasToKeep,:);
    pltData.areaNamesInSamples=pltData.areaNamesInSamples(areasToKeep);
end



if verbose
    fprintf('\n\n')
end
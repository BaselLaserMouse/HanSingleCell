function varargout = listAllTracedCells(inARA)
% return a cell array listing all traced cells
%
% function [allCells,summary,byAnimal,byCell] = listAllTracedCells(inARA)
%
% Purpose
% Return map structure of all cells that have been traced and exported to the ARA space. 
%
%
% Usage:
% cd to experiment root directory and run function listAllTracedCells
% output argument is optional and if present we don't print to screen. 
%
% Inputs
% If inARA is true (true by default) then return paths to cell traced
% in ARA space. Otherwise return paths to cell traces in sample space.
%
% If inARA is a string (unlikely), we look for traced cells in this relative
% path in each sample directory. i.e. when inARA is true what the code 
% does is look for traced cell data in 'downsampled/sample2ARA'. So setting
% inARA to true is the same as setting it to 'downsampled/sample2ARA' 
% 
%
% Outputs
% allCells - a map where keys are cell IDs and values are paths to CSV files of trees.
% summary - summary output of aratools.utils.returnSparseDataSummary
% byAnimal - byAnimal output of aratools.utils.returnSparseDataSummary
% byCell - byCell output of aratools.utils.returnSparseDataSummary
%
%
% Rob Campbell - Basel 2016
%
% also see:
% aratools.utils.returnSparseDataSummary



if nargin<1
    inARA=1;
end

if inARA==1
    subDir=fullfile('downsampled','sample2ARA');
elseif inARA==0
    subDir='downsampled';
elseif isstr(inARA)
    subDir=inARA;
end

[dirs,details] = aratools.utils.returnProcessedExperiments(subDir); %Read csv file "analysis_log.csv" and process all cells listed within it

if isempty(dirs)
    fprintf('%s found no directories with pattern %s\n',mfilename,points)
    return
end
if nargout<1
    fprintf('Found %d data directories containing these cells:\n\n',length(dirs))
end

[summary,byAnimal,byCell]=aratools.utils.returnSparseDataSummary; %Build the cell database

allCells={};
cellID={};
for ii=1:length(dirs)
        csvFiles = dir(fullfile(dirs{ii},'*.csv'));
        if isempty(csvFiles) 
            fprintf('** Found no CSV files in directory %s  SKIPPING\n',dirs{ii})
            continue
        end

        thisAnim = byAnimal.(details(ii).animalID);

        %Build the expected names of the CSV files and try to load them. 
        traces = thisAnim.traces;

        for kk=1:length(traces)     
            fname = sprintf('%s_traced_cells_tree_%02d.csv', thisAnim.animal, traces(kk).cellNum);

            if traces(kk).excludeFromAnalysis
                fprintf('Excluding cell %s from analysis\n', fname)
                continue
            end
            pathToCSV = fullfile(dirs{ii},fname);
            if ~exist(pathToCSV,'file')
                fprintf('Failed to find csv file %s. SKIPPING!\n', pathToCSV)
                continue
            end

            if nargout<1
                fprintf('%s\n',  pathToCSV)
            end
            allCells{length(allCells)+1}=pathToCSV;         
            cellID{length(cellID)+1}=thisAnim.traces(kk).cellID;

        end

end %for ii=1:length(dirs)

if nargout>0
    varargout{1} = containers.Map(cellID,allCells);
end
if nargout>1
    varargout{2}=summary;
end
if nargout>2
    varargout{3}=byAnimal;
end
if nargout>3
    varargout{4}=byCell;
end

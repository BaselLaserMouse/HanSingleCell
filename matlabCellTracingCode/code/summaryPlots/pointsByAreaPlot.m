function varargout=pointsByAreaPlot(thisXylem,varargin)
% Plot matrix of projection pattern of neurons by brain area
%
%
% function pltData=pointsByAreaPlot(thisXylem,'Param1','Value1',...)
%
% Purpose
% This is a critical data analysis plot. It plots the number of points
% by brain area for the population of sample brains or sample neurons. 
% The function receives the output pointsInARA. This can be a vector
% of such structures. If the structure has been processed by getLeavesInARA,
% then the notes field says 'leaves' and the data are treated accordingly 
% (first index is the root node and the rest area leaves).
%
%
% ----------------------------------------------------------------------------------
% Inputs (required)
% thisXylem- a xylem data structure. This is an array of outputs from getLeavesInARA 
%        or pointsInARA that has been transformed into the xylem class. e.g.
%        Xdata = xylem(data)
%
%
%
% ----------------------------------------------------------------------------------
% Inputs (optional)
%  This function has *MANY* optional input arguments. In the interests of sanity 
%  we will group these below by category. 
%
%
% - Arguments that relate to the sort of metric being plotted
% * 'dataMetric'           - ['leaves' by default] Other options are the field names 
%                            in thisXylem.data(1).pointsInARA. If
%                            you want to plot axon length by area
%                            use 'upSampledPoints'. This will return the number of 
%                            mm in each projected area for each cell.
%                             The data field is the output of aratools.sparse.assignToARA. 
% * 'doLog'                - [optional, true/false, false by default] does a log10 of the 
%                             matrix before plotting.
% * 'normaliseBySampleSum' - false by default. if true, we normalise each cell by the sum 
%                            of it's connections.
%
%
% - Arguments relating to the exclusion or inclusion of different brain areas
% * 'showAreas' -  a cell array of strings indicating which areas to plot and in what order.
%                  These are the long names. The string 'Primary visual area' will match
%                  'Primary visual area, layer 4', 'Primary visual area, layer 6', etc...
%                  By default this is missing and all are shown.
%                  If provided, the areas are arranged in the order listed in this cell 
%                  array unless the user sets sortAreasByStrength to true.
%                  ** See the helper function: brainAreaNames.visualAreas **
% * 'excludeAreas' - a cell array of strings indicating which areas to
%                  exclude from the analysis. These are the long 
%                  names. The string 'Primary visual area' will match
%                  'Primary visual area, layer 4', 'Primary visual area, layer 6', etc...
%                  NOTE: this just does a string-based match. No reference
%                  to the hierarchical organization of the tree is made. 
%                  By default this is missing and none are removed.
%                  The excludeAreas operation is always run after showAreas.
%                  e.g. so you can add all of V1 then remove L1.
%
% * 'excludeBorders' - Disabled by default. Otherwise this is a scalar that defines
%                      which points are going to be excluded based on their proximity to a 
%                      a border that is not white matter or out of the brain. e.g. if
%                      excludeBorders is 3 then all points within 3 voxels microns of a border
%                      do not contribute to the plot. 
% * 'excludeSomataNonV1' - true by default. If true, the matrix does not contain cells with
%                          somata outside of V1 layer 2/3. Also excluded are cells that were
%                          labeled with "forceNotV1: 1" in the details file
%
%
% - Arguments that relate to the order of the data in the plot
% * 'sortAreasByStrength'       - [false by default] if true, sort areas by connection strength
% * 'sortSamples'               - [false by default] if true, we do hierarchical clustering and order
%                                 the neurons according to this. Not done if groupSamples is not empty.
%                                 alternatively, this can be the clusterStats field from a different 
%                                 plot's pltData. Then the samples are sorted according to this. 
% * 'exludeAreasFromClustering' - [empty by default] Otherwise a cell array of area names
% 
% 
% - Arguments that are somewhat ad hoc and we likely don't need
% * 'segregateL234' - [optional, true/false, false by default] if true, the function
%                     pushes all neurons with somata not in L2/3 or L4 to columns on the 
%                     far right of the plot. Note the groupSamples takes precedence over this.
% * 'groupSamples'  - a cell array of cell arrays of strings or a cell array of strings. 
%                     If it's a cell array of cell arrays then we demarcate the sub-groups visually.
%                     If present, only these samples are kept. Must be neuron ID names e.g. YH123-01
%                     If present, segregateL234 and sortSamples are set to false.
% * 'includeBothHemispheres' - [false be default] If false, we plot only points on the same hemisphere
%                              as the soma. 
%
% - Arguments loosely relating to plot tidiness 
% * 'suppressYTickLabels' - [false by default]
% * 'suppressSomata'      - [true by default]
% * 'saveName'            - [empty by default] Otherwise if it's a string, it's saved 
%                           to the current dir with this name.
%
%
%
% % ----------------------------------------------------------------------------------
% Important arguments that are part of xylem (optional)
%
% The xylem class is responsible for filtering and pre-processing of the data. 
% By "filtering" we mean operations such as pooling from multiple layers into a single 
% area, excluding points near area boundaries, etc. You can set what is done using
% properties of xylem. Important parameters include:
%
% xylem.groupLayers - boolean
% xylem.groupAreas -  a string or cell array corresponding to valid grouping 
%                actions. See help aratools.sparse.groupARAindexes
%
%
%
%
% ----------------------------------------------------------------------------------
%
% Outputs
% pltData - a structure containing the data used to make the plot. 
% H - plot handles
%
%
%
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
% Examples
%
% ONE - Plot non-logged data with three areas removed:
%  pointsByAreaPlot(data,'doLog',false,'excludeAreas',{'out','corpus cal', 'fiber tracts'})
%
%
% TWO - Plot only a selection of visual areas and do so in the order specified in the cell array, "pltAreas"
%       [pooled are the data grouped by layers - see groupARAindexes.m]
%  pointsByAreaPlot(pooled,'doLog',false,'showAreas',brainAreaNames.visualAreas,'sortSamples',false)
%
%
% THREE - Do not sort columns by hierarchical clustering
%  pointsByAreaPlot(pooled,'sortSamples',false)
%
%
% FOUR - Enforce a particular ordering on the samples (cells/columns)
% Note, can be a single cell arrays of strings
%  group{1}={'YH141-01','YH130-02','YH175-03'};
%  group{2}={'YH124-01','YH134-01'};
%  group{3}={'YH84-01','YH174-01'};
%  group{4}={'YH140-01','YH175-02','YH169-01'};
%  group{5}={'YH163-01','YH176-01','YH176-02','YH167-01'};
%
%  pointsByAreaPlot(pooled,'groupSamples',group)
% 
%
%
% Rob Campbell - Basel 2015
%
%
% See also: pointsInARA, brainAreaBarChart, groupARAindexes, xylem


if nargin==0
    help(mfilename)
    return
end

if ~isa(thisXylem,'xylem')
    fprintf('input "thisXylem" must be of class xylem.\nOnce you have your processed tree data you should do, for instance:\n x=xylem(data)\n')
    return
end



% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
%Handle input arguments
params = inputParser;
params.CaseSensitive = false;
params.addParamValue('dataMetric','leaves',@ischar)
params.addParamValue('doLog',false,@islogical)
params.addParamValue('segregateL234',false,@islogical)
params.addParamValue('sortAreasByStrength',false,@islogical)
params.addParamValue('sortSamples',false,@(x) islogical(x) || isstruct(x))
params.addParamValue('excludeAreas',{},@iscell) 
params.addParamValue('excludeBorders',-1,@isscalar) 
params.addParamValue('groupSamples',{},@iscell)
params.addParamValue('showAreas',{},@iscell)
params.addParamValue('excludeAreasFromClustering',{},@iscell)
params.addParamValue('normaliseBySampleSum',false,@islogical)
params.addParamValue('suppressYTickLabels',false,@islogical)
params.addParamValue('suppressSomata',true,@islogical)
params.addParamValue('saveName','',@ischar)
params.addParamValue('includeBothHemispheres',false,@islogical)
params.addParamValue('excludeSomataNonV1',true,@islogical)

params.parse(varargin{:});

doLog=params.Results.doLog;
dataMetric=params.Results.dataMetric;
excludeAreas=params.Results.excludeAreas;
showAreas=params.Results.showAreas;
groupSamples=params.Results.groupSamples;
segregateL234=params.Results.segregateL234;
sortAreasByStrength=params.Results.sortAreasByStrength;
sortSamples=params.Results.sortSamples;
excludeAreasFromClustering=params.Results.excludeAreasFromClustering;
excludeBorders=params.Results.excludeBorders;
normaliseBySampleSum=params.Results.normaliseBySampleSum;
saveName=params.Results.saveName;
suppressYTickLabels=params.Results.suppressYTickLabels;
suppressSomata=params.Results.suppressSomata;
includeBothHemispheres=params.Results.includeBothHemispheres;
excludeSomataNonV1=params.Results.excludeSomataNonV1;


if isstruct(sortSamples)
    clusterStats=sortSamples;
    sortSamples=true;
else
    clusterStats=[];
end



if ~sortSamples
    if segregateL234
        fprintf('sortSamples is false. Setting segregateL234 to false\n')
    end
    segregateL234=false;
    sortSamples=false;
end


% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 




%Get lists of area IDs and names, etc
%if we have leaves, then the first leaf in each case is the location of the cell body

if strcmpi(dataMetric,'leaves') || strcmpi(dataMetric,'upSampledPoints')
    startInd=2; %Gets rid of the soma point
else
    startInd=1;
end
pltData = GetPointsByArea(thisXylem,startInd,dataMetric,excludeBorders,includeBothHemispheres,excludeSomataNonV1);

if startInd==2
    somaLocationInPlot = zeros(length(pltData.dataStructure),1);
    localIncomplete = zeros(length(pltData.dataStructure),1);
    for ii=1:length(pltData.dataStructure)
        somaLocationInPlot(ii) = pltData.dataStructure(ii).pointsInARA.(dataMetric).ARAindex(1);
        localIncomplete(ii) = pltData.dataStructure(ii).details.localIncomplete;
    end
else
    localIncomplete=[];
    somaLocationInPlot=[];
end




%If the user asked for a subset of the rows (areas) only, then we perform this
%operation here
if ~isempty(showAreas)
    fprintf('Keeping selected brain areas\n')
    pltData = filterPlotAreasFromPltData(pltData,showAreas,1);

    f=find(sum(pltData.dataMat)==0);
    if ~isempty(f)
        fprintf('WARNING: you have %d cells with no hits in the selected brain areas\n', length(f))
    end
end

%Exclude if needed
if ~isempty(excludeAreas)
    fprintf('Peforming area exclusion\n')
    pltData = filterPlotAreasFromPltData(pltData,excludeAreas);
end


pltData=sortAreas(pltData,sortAreasByStrength); %sort area names

%Add these variables here (after sortAreas)
pltData.localIncomplete=localIncomplete;
pltData.somaLocationInPlot = somaLocationInPlot;


%sort samples by their hierarchical clusters
if length(pltData.cellIDs)>5 & sortSamples
    
    if 0 % permute labels for each cell separately PERMUTATION TEST (TODO: to properly elsewhere)
        for ii=1:size(pltData.dataMat,1)
            R=randperm(size(pltData.dataMat,2));
            pltData.dataMat(ii,:) = pltData.dataMat(ii,R) ;
        end
    end

    if isempty(clusterStats)
        fprintf('--> Sorting samples by similarity using hierarchical clustering\n')
        clusterStats=clusterTracedCells(pltData,excludeAreasFromClustering);
    else
        fprintf('Reordering using a previously calculated clustering run\n')
    end

    pltData = reorderCols(pltData,clusterStats.perm);
    pltData.clusterStats=clusterStats;
end



%Sort samples
if ~isempty(groupSamples)
    fprintf('Sorting samples\n')
    if iscell(groupSamples{1})
        samplesToKeep = [groupSamples{:}];
    else
        samplesToKeep = groupSamples;
    end

    %re-order based on the above
    ind=[];

    for ii=1:length(samplesToKeep)
        thisS=samplesToKeep{ii};

        thisInd = strmatch(thisS,pltData.cellIDs);
        if isempty(thisInd)
            fprintf(' * Unable to find %s in dataset. Skipping! * \n',thisS)
            continue
        end

        ind=[ind,thisInd];

    end

    %re-build the data with these groups
    pltData = reorderCols(pltData,ind);

end




if segregateL234 
    fprintf('Segregating Layer 234')
    structures = getAllenStructureList; %get the Allen structure names
    f23=strmatch('Primary visual area, layer 2/3',structures.name);
    f4=strmatch('Primary visual area, layer 4',structures.name);
    
    L23=structures.id(f23);
    L4=structures.id(f4);

    f=(pltData.somaLocation~=L23 & pltData.somaLocation~=L4);
    
    ind=1:length(pltData.somaLocation);
    ind = [ind(~f),ind(f)];

    pltData = reorderCols(pltData,ind);
end



%--------------------------------------------
%plot
clf


if doLog
    fprintf('Doing log10(dataMat)+0.5\n')
    pltData.dataMat=log10(pltData.dataMat+0.5);
end


im=pltData.dataMat'; %TODO: this rotation of the plot matrix is potentially problematic. Is this the best way to go here?
if normaliseBySampleSum
    sumOfRows = sum(im,2);
    im=im ./ repmat(sumOfRows,1,size(im,2));
end
imagesc(im)



%Set up the color map
if ~doLog
    r=range(pltData.dataMat(:))*100; %Make sure we have more than enough values in the scale or very low values aren't plotted
    colormap(parula(r)) %Important: without specifying the range we loose the 1s
elseif doLog
    r=length(unique(pltData.dataMat(:)));
    colormap(parula(r))
end
%Make zero connections 4/5 black
C=colormap;
C(1,:) = [0.25,0.25,0.25];
colormap(C)


%add soma location for leaf data
if strcmpi(dataMetric,'leaves')
    hold on
    for ii=1:length(pltData.somaLocationInPlot)

        f=find(pltData.allAreas(:,1)==pltData.somaLocationInPlot(ii));
        if isempty(f)
            if ~suppressSomata %warn the user if we're missing somata areas from the plot
                fprintf('You seem to have removed the area containing the cell body of %s.\n',pltData.cellIDs{ii})
            end
            continue
        end

        H.soma(ii)=plot(f,ii,'r*', 'MarkerSize', 8, 'LineWidth', 1.5);
        if pltData.localIncomplete(ii)
            set(H.soma(ii),'Color',[1,0.5,0.5])
        end

    end
    hold off
end


%Let's check whether we have duplicate area names and issue a warning if this is taking place
aNames = pltData.areaNamesInSamples;
[U,~,idx] = unique(aNames);
if length(U) ~= length(aNames)
    fprintf('There are duplicate area names!\n')
    unique_idx = accumarray(idx(:), (1:length(idx))',[] , @(x) {sort(x)}); %applies function each subset of second input arg

    for ii=1:length(unique_idx)
        if length(unique_idx{ii})>1
            fprintf('%s\n', aNames{unique_idx{ii}(1)})
        end
    end

end


%Overlay start and end points of area groups
b=pltData.groupingPoints;

if ~isempty(b)
    hold on
    for ii=1:length(b)
        H.groupingLinesW(ii)=plot(xlim,[b(ii),b(ii)]+0.5,'w--');

        H.groupingLinesK(ii)=plot([0.5,-1.2],[b(ii),b(ii)]+0.5,'k--');
    end
    hold off
    set(H.groupingLinesK,'Clipping','off')
end


%lines between columns
hold on
for ii=1:length(pltData.areaNamesInSamples)-1
    H.colSepLines(ii)=plot([ii,ii]+0.5,ylim,'-k','LineWidth',1);

end


%Overlay dashed white lines if we have sub-groups defined
if ~isempty(groupSamples) & iscell(groupSamples{1})
    s=cellfun(@length,groupSamples);
    if sum(s) == size(pltData.dataMat,2) %otherwise we're missing cells and the boundaries will be wrong
        for ii=1:length(s)-1
            x=sum(s(1:ii))+0.5;
            y=ylim;
            H.sampleGroupingLines(ii)=plot([x,x],y,'w--');
        end
    end

end
hold off


%Add a colorbar
H.colorbar=colorbar;
if doLog
    colString = 'number of connections (colors logged)';
else
    colString = 'number of connections';
end



%Arrange the labels 
set(get(H.colorbar,'XLabel'),'String',colString)
if doLog
    xtick = get(c,'XTick');
    set(H.colorbar,'XTick',[0,xtick(end)/2,xtick(end)])
end

set(gca,'XTick',1:size(pltData.dataMat,1), 'XTickLabel',pltData.areaNamesInSamples,...
    'YTick',1:size(pltData.dataMat,2), 'YTickLabel',regexprep(pltData.cellIDs,'_','-') )

if suppressYTickLabels
    set(gca,'YTick',[])
end

if verLessThan('matlab','8.4')
    fprintf('You are running <R2014b. Can not rotate tick labels.\n')
else
    set(gca,'XTickLabelRotation',45)
end

if suppressSomata && isfield(H,'soma')
    delete(H.soma)
end





%TODO: decide if saving plot should be done here and if so how best to do it
if ~isempty(saveName)
    fsize=10;
    set(gcf,'PaperPosition',[0,0,10,10])
    set(gca,'FontSize',fsize)
    print('-depsc',[saveName,'.eps'])
    print([saveName,'.png'],'-dpng','-r200')
end



% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
if nargout>0
    pltData.args = params.Results;
    varargout{1}=pltData;
end

if nargout>1
    varargout{2}=H;
end




function pltData = reorderCols(pltData,perm)
    fprintf('Reordering cells\n')
    pltData.dataMat = pltData.dataMat(:,perm);
    pltData.cellIDs = pltData.cellIDs(perm);
    pltData.localIncomplete = pltData.localIncomplete(perm);
    pltData.somaLocationInPlot = pltData.somaLocationInPlot(perm);
    pltData.somaLocation = pltData.somaLocation(perm);



function pltData=GetPointsByArea(thisXylem,startInd,structField,excludeBorders,includeBothHemispheres,excludeSomataNonV1)
    % This is function returns points or leaves by brain area
    % 

    data = thisXylem.returnData('excludeBorders',excludeBorders);


    %Remove cells which don't have a soma in V1 L2/3
    if excludeSomataNonV1
        SS=scc.diagnostic.summariseTerminationTypes(thisXylem,true);
        nonV1CellIDs = SS(5).IDsOfAllCells;
        fprintf('\nRemoving the following %d neurons which don''t have somata in V1 layer 2/3\n',length(nonV1CellIDs))
        for ii=length(data):-1:1
            if ~isempty(strmatch(data(ii).details.cellID,nonV1CellIDs))
                fprintf('Removing %s\n', data(ii).details.cellID)
                data(ii)=[];
            end
        end
    end



    IDs = []; %vector containing a list of area IDs that have leaves or nodes
    for ii=1:length(data)
        IDs = [IDs;data(ii).pointsInARA.(structField).ARAindex(startInd:end)];
    end


    allAreas = unique(IDs); %These are the only areas present in these data sets

    fprintf('%s found %d samples that together contain points in %d unique brain areas.\n', ...
             mfilename, length(data), length(allAreas))


    verbose=1;


    %Count the number of hits per area
    allAreas(:,2) = 0;
    for ii=1:size(allAreas,1)
        f = find(IDs(:,1)==allAreas(ii,1));
        allAreas(ii,2) = length(f);
    end

    %Get the names of the brain areas
    areaNamesInSamples = cell(size(allAreas));
    for ii=1:length(allAreas)
        areaNamesInSamples{ii} = structureID2name(allAreas(ii,1));
    end

    %Make the data matrix which we will plot
    dataMat = zeros(size(allAreas,1),length(data));

    verbose=0;
    if verbose
        fprintf('\n\nStarting data extraction in GetPointsByArea\n\n')
    end

    cellIDs={};
    for ii=1:length(data)
        numHits=0;

        if verbose
            fprintf(' --> Sample %s\n',data(ii).details.cellID)
        end

        cellIDs{ii} = data(ii).details.cellID;
        for jj=1:length(allAreas)


            if includeBothHemispheres==true
                ARAindex = data(ii).pointsInARA.(structField).ARAindex;
            elseif includeBothHemispheres==false
                ARAindex = data(ii).pointsInARA.(structField).ARAindex;
                hemisphere = data(ii).pointsInARA.(structField).hemisphere;
                f=find(hemisphere==hemisphere(1));
                ARAindex=ARAindex(f);
            end
                

            f = find( ARAindex == allAreas(jj) );
            numHits = numHits+length(f);
            if isempty(f), continue, end
            dataMat(jj,ii) = length(f);

            if verbose
                fprintf('%d hits in area %d: %s\n',length(f),allAreas(jj), areaNamesInSamples{jj,1})
            end 

        end 
        if verbose
            fprintf('\nFound %d hits in %s\n ', numHits, data(ii).pointsFname)
        end
    end

    %If we're using up-sampled points we should convert this value in voxels to mm
    if strcmp(structField,'upSampledPoints')
        %Get the number of microns per pixel from the number of voxels per pixel
        upRes = data(1).pointsInARA.upSampledPoints.details.upSampleResolution; %Number of voxels per pixel
        micronsPerVoxel = str2double(data(1).voxelSize);
        micronsPerUpSampledPoint = upRes * micronsPerVoxel;
        dataMat = dataMat * micronsPerUpSampledPoint;
        dataMat = dataMat * 1E-3; %Convert to mm
    end


    %These will be used as outputs and for plotting
    pltData.dataMat = dataMat;
    pltData.areaNamesInSamples=areaNamesInSamples(:,1);
    pltData.allAreas = allAreas;
    pltData.dataStructure = data;
    pltData.cellIDs=cellIDs;
    pltData.somaLocation = [data.somaIndex];
function varargout = getLaminarData(cleanCells,areaName,atlas,dataMetric)
    % Obtain distributions of terminals by layer within a subset of defined areas
    %
    % out = getLaminarData(cleanCells,areaName,atlas,dataMetric)
    %
    % Find the distribution of axon in a named target area and return as a structure.
    %
    % Inputs
    % cleanCells - Xylem data structure containing the selected good cells (e.g. no crappy 
    %              ones where the brain was falling apart)
    % areaName - Name of the area to plot (full string e.g. "Primary visual area")
    % atlas - output of aratools.atlascacher.getCachedAtlas (optional, done automatically 
    %         by default but takes ~2.5 s per function call)
    %
    % Outputs
    % out - a structure that contains the data that can be plotted and details of what the 
    %       data contain.
    %
    % Examples
    % >> load ~/tvtoucan/Mrsic-Flogel/hanyu/Analyses/cleanCells.mat
    % >> d=getLaminarData(cleanCells,'Primary visual cortex');
    %
    % >> visNames=brainAreaNames.visualAreas;
    % >> clear OUT
    % parfor ii=1:length(visNames); disp(ii), OUT(ii)=getLaminarData(cleanCells,visNames{ii}, LOADED_ARA); end
    %
    %
    % Rob Campbell - Basel 2017


    if nargin<2 || isempty(areaName)
        areaName = 'Lateral visual area';
    end

    if nargin<3 || isempty(atlas)
        atlas = aratools.atlascacher.getCachedAtlas;
    end

    if nargin<4 
        dataMetric = 'upSampledPoints';
    end





    %Find which cells project where and remove all projections with under 1 mm of axon
    D=pointsByAreaPlot(cleanCells,'dataMetric', 'upSampledPoints', 'excludeBorders', 2, 'excludeSomataNonV1', true);
    M = D.dataMat;
    M(M<1)=0; % Exclude areas with less than 1 mm of axon

    %If user asked for leaves we re-do and filter the matrix by M, above
    if strcmp(dataMetric,'leaves')
        D=pointsByAreaPlot(cleanCells,'dataMetric', 'leaves', 'excludeBorders', 2, 'excludeSomataNonV1', true);
        M(M>0)=0; %We will filter by this
        M = D.dataMat .* M;
    end


    % Find cells projecting to the area in question 
    ind=strmatch(areaName,D.areaNamesInSamples);

    if isempty(ind)
        fprintf('Can not find area %s in projection list\n', areaName)
        return
    end

    thisArea = find(M(ind,:)>0);


    [childAreas,parentInd] = getChildAreas(areaName);
    data = cleanCells.returnData('excludeBorders',2);
    details = [data.details];
    cellIDs = {details.cellID};
    thisCellID={};

    COUNTS=[];
    PLANES=[];
    traceFiles = {}; %The location of the raw data file
    for ii=1:length(thisArea)
        thisCell = D.cellIDs{thisArea(ii)};
        thisInd = strmatch(thisCell,cellIDs);
        theseData = data(thisInd);
        traceFiles{ii} = theseData.pointsFname;
        thisCellID{ii} = thisCell;
        [tmpC,tmpP]=countsInLayers(theseData,childAreas,parentInd,atlas);
        COUNTS=cat(3,COUNTS,tmpC); %The number of points in the up-sampled axons
        PLANES = [PLANES,tmpP]; %The planes in the atlas that contain axon in this cortical are
    end

    if nargout>0
        out.counts = COUNTS;
        out.areaTable = childAreas;
        out.areaName = areaName;
        out.cellID = thisCellID;
        out.planesInARA = unique(PLANES);
        out.traceFiles=traceFiles;
        out.dataMetric=dataMetric;
        varargout{1}=out;
    end




% Internal functions follow 
function [out,planesWithAxonInThisArea]=countsInLayers(theseData,childAreas,parentInd,atlas)
    % out - n by 2 with first column being ID and second being counts

    f = find(theseData.pointsInARA.upSampledPoints.ARAindex==parentInd);
    thesePoints = theseData.pointsInARA.upSampledPoints.ARAindex(f,:);


    %get the index of each of these
    inds = ones(1,length(f));

    if length(inds)==0
        fprintf('Failed to find area %d in cell %s\n', parentInd, theseData.cellID)
        out=[];
        return
    end

    planesWithAxonInThisArea=[];
    for ii=1:length(f)
        sp = theseData.pointsInARA.upSampledPoints.sparsePointMatrix(f(ii),:);
        inds(ii) = atlas.atlasVolume(sp(1),sp(2),sp(3));
        planesWithAxonInThisArea(ii)=sp(3);
    end

    [counts,IDs] = hist(inds, unique(inds));

    % Order these neatly so that each layer has a count value attached to it
    allLayerIDs = childAreas.id;
    allCounts = zeros(size(allLayerIDs));

    % Add up the number of voxels in this layer including only planes where the
    % neuron has axon in this area. 
    thesePlanes = atlas.atlasVolume(:,:,planesWithAxonInThisArea);
    for ii=1:length(allLayerIDs)

        voxelsInLayer(ii) = length(find(thesePlanes==allLayerIDs(ii)));

        f=find(IDs==allLayerIDs(ii));
        if isempty(f), continue, end
        allCounts(ii) = counts(f);

    end
    
    out = [allLayerIDs,voxelsInLayer',allCounts];
    planesWithAxonInThisArea = unique(planesWithAxonInThisArea);

function [childTable,parentInd] = getChildAreas(areaName)
    childTable = getAllenStructureList('childrenOf',areaName);
    % remove parent area
    ind = strmatch(areaName,childTable.name,'exact');
    if ~isempty(ind)
        parentInd = childTable.id(ind);
        childTable(ind,:)=[];
    end

    [~,ind]=sort(childTable.name);
    childTable=childTable(ind,:);

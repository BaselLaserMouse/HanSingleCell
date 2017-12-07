function varargout = clusterPos(cleanCells,highlightName)
    % Highlight soma positions of cells projecting to a given brain area
    %
    % [plotHandles, dataTable] = clusterPos(cleanCells,highlightName)
    %
    % Purpose
    % Plot the locations of somata that project to a given area
    %
    % Inputs
    % cleanCells - The "clean cells" as a Xylem data object.
    % highlightName - index of visual area to which to plot projections (a number from 1 to 13)
    %
    %
    % e.g.
    % >> load ~/tvtoucan/Mrsic-Flogel/hanSingleCell2017/Analyses/cleanCells.mat
    % >> clusterPos(cleanCells,3)
    %
    %
    % % Go through all projection targets and make an output table for each
    % >> for ii=2:13, [~,dt{ii}]=clusterPos(cleanCells,ii); end
    %
    %
    % Rob Campbell - Basel 2017


    if nargin<2
        highlightName=[];
    end


    A=aratools.atlascacher.getCachedAtlas;
    load('brainAreaProjections.mat') %From: ReferenceAtlas/ARA_CCFv3/ARA_25_micron_mhd/projections/

    %Make the projection matrix data
    D=pointsByAreaPlot(cleanCells,'dataMetric', 'upSampledPoints', 'excludeBorders', 2, 'excludeSomataNonV1', true);
    D.dataMat(D.dataMat<1)=0; %Exclude areas with less than 1 mm of axon

    %The area names and abbreviations:
    [n,c,abrv]=brainAreaNames.visualAreas; 

    if ~isempty(highlightName) && isnumeric(highlightName)
        highlightName = c.areaNames{highlightName};
        fprintf('Highlighting cells projecting to %s\n', highlightName);
    end

    %Generate plot that we will modify
    H=overlayCellsOnThreeProjections([],projections,c);


    delete(H.axesCoronal)
    H = rmfield(H,'axesCoronal');
    H = rmfield(H,'coronal');

    delete(H.axesSagittal)
    H = rmfield(H,'axesSagittal');
    H = rmfield(H,'sagittal');

    set(H.axesTransverse,'position',[0.05,0.05,0.8,0.9])


    % Make all highlighted areas areas gray apart from the selected one
    if ~isempty(highlightName)
        for ii=1:length(H.transverse.areas)
            if strcmpi(H.transverse.areas(ii).Tag, 'AREA BOUNDARY') || ...
                strcmpi(H.transverse.areas(ii).Tag, 'Primary visual area') || ...
                strcmpi(H.transverse.areas(ii).Tag,  highlightName)
                continue
            end 

            set(H.transverse.areas(ii),'EdgeColor',[1,1,1]*0.25, 'FaceColor',[1,1,1]*0.75)
        end
    end

    % labels
    hold on

    params={'FontSize', 18, 'FontWeight', 'Bold'};

    H.tRS=text(276,180,'RS (lateral)');
    set(H.tRS,'rotation',90, params{:})

    %H.tV1=text(330,165,'V1');
    %set(H.tV1,'rotation',0, params{:})

    H.tA1=text(420,210,'A1');
    set(H.tA1,'rotation',0, params{:})

    H.tS1=text(367,236,'S1');
    set(H.tS1,'rotation',0, params{:})

    H.tPM=text(292,180,'PM');
    set(H.tPM,'rotation',63, params{:})

    H.tAM=text(291,217,'AM');
    set(H.tAM,'rotation',17, params{:})

    H.tA=text(323,229,'A');
    set(H.tA,'rotation',0, params{:})

    H.tRL=text(354,221,'RL');
    set(H.tRL,'rotation',-33, params{:})

    H.tAL=text(379,197,'AL');
    set(H.tAL,'rotation',0, params{:})

    H.tLM=text(383,168,'LM');
    set(H.tLM,'rotation',-90, params{:})

    H.tLI=text(396,168,'LI');
    set(H.tLI,'rotation',-90, params{:})

    H.tP=text(373,124,'P');
    set(H.tP,'rotation',0, params{:})

    H.tPOR=text(393,132,'POR');
    set(H.tPOR,'rotation',35, params{:})

    H.tTEA=text(411,159,'TEA');
    set(H.tTEA,'rotation',35, params{:})



    % Plot soma positions
    if ~isempty(highlightName)
        set(gcf,'Name',(sprintf('Highlighting cells projecting to "%s"\n', highlightName)) );
    end
    xlim([260,440])
    ylim([100,250]) 


    %Overlay the bisection lines 

    r_h=111; %circle radius
    c_x=442; %circle offset in x
    c_y=177; %circle offset in y


    x=linspace(0,1,50)+2.75; % X positions to take the sin and cos of for plotting
    V1divmrkr = {'--','LineWidth',1.5,'Color',[1,1,1]*0.2};
    plot(c_x+cos(x)*r_h, c_y+sin(x)*r_h, V1divmrkr{:})

    %The coefs for the line
    b0_v=87.0;
    b1_v=0.2353;
    x=[289,374];
    plot(x, b0_v + x*b1_v, V1divmrkr{:})


    %Remove non-V1 cells
    S=scc.diagnostic.summariseTerminationTypes(cleanCells,true);


    %Plot the somata
    dataFromCleanCells = cleanCells.returnData;


    mSize=12;

    %Plot non-abrupt terminations of V1 non-back-labelled cells
    cleanCellIDs = S(4).IDsofCleanCells;
    [cellLocs(1),cm(1)]=plotSomaPositions(cleanCellIDs,dataFromCleanCells,D,highlightName,{'ob','MarkerSize',mSize},...
        {'ob', 'MarkerFaceColor', [0.5,0.5,1],'MarkerSize',mSize});

    % Plot cells with abrupt terminations (V1 non-back-labelled cells)
    details = [dataFromCleanCells.details];
    details=details(~S(4).cleanCell & S(4).allCells);
    prematureCellIDs = {details.cellID};
    [cellLocs(2),cm(2)]=plotSomaPositions(prematureCellIDs,dataFromCleanCells,D,highlightName,{'sb','MarkerSize',mSize+1},...
          {'sb', 'MarkerFaceColor', [0.5,0.5,1],'MarkerSize',mSize+1});


    % Now we make a table that contains the following columns

    % 1. Cell ID (string)
    % 2. Projection target name (string) - This will be the target for this run of clusterPos unless it's called with no arguments (in which case all are returned)
    % 3. Medial or Lateral ('M' or 'L') side of the border
    % 4. Rostral or caudal ('R' or 'C') side of the border
    % 5. Is premature? [bool]
    % 6. Soma position

    IDs = cat(2,cm.id)';
    L = length(IDs);

    %Figure out on which side of each boundary each cell falls;
    %Cells within the circle diamater will be lateral so make all medial by default
    ML = repmat({'M'},L,1);

    pos = cat(1,cellLocs.pos);

    d = squareform(pdist([c_x,c_y;pos])); %Distance of all WRT circle centre
    d = d(1,2:end);
    lateralInds = find(d<=r_h);
    ML(lateralInds) = repmat({'L'},1,length(lateralInds));


    %Cells beneath the line will be caudal. Make this the default then figure out what is rostral.
    RC = repmat({'C'},L,1);
    border=b0_v + pos(:,1)*b1_v;
    rostralInds = find((pos(:,2)-border)>0);
    RC(rostralInds) = repmat({'R'},1,length(rostralInds));

    %Cells that project to the area in question
    proj = repmat({''},L,1);

    if ~isempty(highlightName)
        projCells = cat(1,cellLocs.projectingCellID);
        for ii=1:length(projCells)
            f=find(strcmp(projCells{ii},IDs));
            proj{f} = highlightName;
        end
    end

    out = table(IDs, ...
        proj, ...
        ML, ...
        RC, ...
        [repmat(false,length(cm(1).id),1); repmat(true,length(cm(2).id),1)], ...
        pos);
    out.Properties.VariableNames = {'CellID','Target','ML','RC','isPremature','position'};

    if nargout>0
        varargout{1} = H;
    end

    if nargout>1
        varargout{2} = out;
    end

    hold off


function [somaLocations,cellMat] = plotSomaPositions(cellIDsToPlot,cellData,projStats,highlightName,allCellsMarker,projectionCellsMarker)
    % return cellMat because some cells might have been removed from it. 

    details=[cellData.details];

    cellID = {details.cellID};

    %Make an output array of cells that have been plotted
    cellMat.id = cellIDsToPlot;

    % Find the index for each cell from the manuscript
    for ii=1:length(cellMat.id)
        f = strmatch(cellMat.id{ii},cellID);

        if ~strcmp(cellData(f).details.cellID, cellMat.id{ii})
            fprintf('WARNING -- clusterPos>plotSomaPosition -- IDs do not match. Something is wrong.\n')
            continue
        end
        % The soma position of each cell in the list in 25 micron voxel units
        cellMat.somaPos{ii} = cellData(f).pointsInARA.rawSparseData.sparsePointMatrix(1,:);
    end

    % Plot all cells
    pltDat = reshape([cellMat.somaPos{:}],3,[])';
    plot(pltDat(:,2), pltDat(:,3),allCellsMarker{:})

    somaLocations.allCells = [(1:length(pltDat))' , pltDat(:,2:3)];
    somaLocations.projectingCells = [];
    somaLocations.pos = pltDat(:,2:3);
    somaLocations.projectingCellID ={};


    % Find cells that project to the highlighted area and overlay them
    if ~isempty(highlightName)
        Hind = strmatch(highlightName,projStats.areaNamesInSamples);

        if isempty(Hind)
            fprintf('Can not find area "%s" in projStats.areaNamesInSamples\n', highlightName)
        else
            Hcells = find(projStats.dataMat(Hind,:)>0); % Index of cells with projections to this area
            HcellIDs = projStats.cellIDs(Hcells); % IDs of cells with projections to this area

            for ii=1:length(HcellIDs)
                tCell = HcellIDs{ii};
                f = strmatch(tCell,cellMat.id);

                if isempty(f) % We will skip some because we cellIDsToPlot will be either premature or non-premature
                    continue
                end

                tmp = cellMat.somaPos{f};
                plot(tmp(:,2), tmp(:,3), projectionCellsMarker{:})

                somaLocations.projectingCells = [somaLocations.projectingCells; tmp];
                somaLocations.projectingCellID =[somaLocations.projectingCellID; tCell];
            end

        end
    end

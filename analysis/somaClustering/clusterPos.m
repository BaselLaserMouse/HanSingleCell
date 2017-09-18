function varargout = clusterPos(cellMat,cleanCells,highlightName)
    % function clusterPos(cellMat,cleanCells,highlightName)
    %
    % Purpose
    % Plot the locations of somata that project to a given area
    %
    % e.g.
    % >> load ~/tvtoucan/Mrsic-Flogel/hanyu/Analyses/cleanCells.mat
    % % The following loads "cellMat"
    % >> load goodCellsMatrixFromJustus.mat   
    % >> clusterPos(cellMat,cleanCells,3)

    if nargin<3
        highlightName=[];
    end


    A=aratools.atlascacher.getCachedAtlas;  
    [visualAreaNames,colorMap]=brainAreaNames.visualAreas;
    load('brainAreaProjections.mat') %From: ReferenceAtlas/ARA_CCFv3/ARA_25_micron_mhd/projections/


    %The area names and abbreviations:
    [n,c]=brainAreaNames.visualAreas; 

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

    H.tV1=text(330,165,'V1');
    set(H.tV1,'rotation',0, params{:})

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

    d = cleanCells.returnData;
    details=cat(1,d.details);
    cellID = {details.cellID};

    % Find the index for each cell from the manuscript
    for ii=1:length(cellMat.id)
        f = strmatch(cellMat.id{ii},cellID);
        cellMat.ind(ii) = f;

        % Soma location in 25 micron voxel units
        cellMat.somaPos{ii} = d(f).pointsInARA.rawSparseData.sparsePointMatrix(1,:);
    end
    pltDat = reshape([cellMat.somaPos{:}],3,[])';

    plot(pltDat(:,2), pltDat(:,3), 'ob')


    %Find cells projected to the highlighted area
    if ~isempty(highlightName)
        Hind = strmatch(highlightName,cellMat.lab);
        if isempty(Hind)
            fprintf('\nERROR: can not find area "%s" in cellmat.lab\n', highlightName)
        end
        Hcells = find(cellMat.mat(:,Hind));
    end

    plot(pltDat(Hcells,2), pltDat(Hcells,3), 'ob', 'MarkerFaceColor', [0.5,0.5,1])


    hold off

    title(sprintf('Highlighting cells projecting to "%s"\n', highlightName));
    xlim([260,440])
    ylim([100,250]) 




    if nargout>0
        varargout{1} = H;
    end

    if nargout>1
        varargout{2} = cellMat;
    end

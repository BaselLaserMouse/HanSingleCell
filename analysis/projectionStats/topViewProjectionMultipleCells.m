function varargout=topViewProjectionMultipleCells(X,laminarData)
    % function varargout=topViewProjectionMultipleCells(cleanCells,laminarDataInstance)
    %
    % Make hairy plot of all defined neurons projecting to one area
    % This is defined by one element of structure, laminarData, which
    % is produced by "getLaminarData"

    clf
    highlightName = laminarData.areaName;

    [~,c]=brainAreaNames.visualAreas;
    load brainAreaProjections

    H=overlayCellsOnThreeProjections([],projections,c);

    delete(H.axesCoronal)
    H.axesTransverse.Position=[0,0.5,0.5,0.5];
    H.axesSagittal.Position=[0,0,0.5,0.5];


    set(H.axesTransverse, 'XLim', [230,435], 'YLim', [105,255])
    set(H.axesSagittal, 'XLim', [105,255],'YLim',[8,180]) 
    hold(H.axesSagittal,'on')
    hold(H.axesTransverse,'on')

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
        for ii=1:length(H.sagittal.areas)
            if strcmpi(H.sagittal.areas(ii).Tag, 'AREA BOUNDARY') || ...
                strcmpi(H.sagittal.areas(ii).Tag, 'Primary visual area') || ...
                strcmpi(H.sagittal.areas(ii).Tag,  highlightName)
                continue
            end 

            set(H.sagittal.areas(ii),'EdgeColor',[1,1,1]*0.25, 'FaceColor',[1,1,1]*0.75)
        end
    end

    % labels


    set(gcf,'Renderer','zbuffer')
    
    hold on 
    H.tracesSagittal={};
    H.tracesTransverse={};

    d = X.returnData;
    details=[d.details];
    cellIDs = {details.cellID};
    drawnow

    %Plot neurite trees
    for ii=1:1:length(laminarData.cellID)
        tCell = laminarData.cellID{ii};
        ind=strmatch(tCell,cellIDs);

        H.tracesTransverse{end+1}=plotTree(d(ind).neuriteTree, H.axesTransverse);
        H.tracesSagittal{end+1}=plotTree(d(ind).neuriteTree, H.axesSagittal);
        if mod(ii,8)==0
            drawnow
        end
    end

    % Overlay somata
    H.somaTransverse={};
    H.somaSagittal={};

    for ii=1:1:length(laminarData.cellID)
        tCell = laminarData.cellID{ii};
        ind=strmatch(tCell,cellIDs);

        H.somaTransverse{end+1}=plotSoma(d(ind).neuriteTree, H.axesTransverse);
        H.somaSagittal{end+1}=plotSoma(d(ind).neuriteTree, H.axesSagittal);
        if mod(ii,8)==0
            drawnow
        end
    end



    set(gcf,'PaperPosition',[0,0,20,20])

    if nargout>0
        varargout{1}=H;
    end



%-------------------------------------------------------------------------------------------
function h=plotTree(data,ax)
    segments = data.getsegments;

    for ii=1:length(segments)
        theseNodes = data.Node(segments{ii});
        theseData = ones(length(theseNodes),3);
        for jj=1:length(theseNodes)
            theseData(jj,:) = theseNodes{jj};
        end

        lineProps={'-','linewidth',0.5, 'Color', [0,0,0.15,0.3]};
        switch ax.Tag
            case 'transverse'
                h(ii)=plot(ax,theseData(:,2), theseData(:,1), lineProps{:});
            case 'sagittal'
                h(ii)=plot(ax,theseData(:,1), theseData(:,3), lineProps{:});
            otherwise
                error('Unknown axes: "%s"\n',ax.Tag)
        end
    end


function h=plotSoma(data,ax)
    %The soma
    theseData = data.Node{1};
    switch ax.Tag
        case 'transverse'
            h=plot(ax,theseData(1,2), theseData(1,1), 'o','MarkerEdgeColor','w','MarkerFaceColor','r');
        case 'sagittal'
            h=plot(ax,theseData(1,1), theseData(1,3), 'o','MarkerEdgeColor','w','MarkerFaceColor','r');
    end


function S=overlaySubCortical(subCorticalInd,ax)
    %
    %
    % Overlay sub-cortical areas on to a view:
    % subCorticalInd - brain area index values 
    %   e.g. for striatum, amygdala and ACC (multiple layers grouped) : {672,131,[211,935,1015]})

    S=[];
    if isempty(subCorticalInd)
        return
    end

    A=aratools.atlascacher.getCachedAtlas;
    vol = A.atlasVolume;
    switch ax.Tag
    case 'transverse'
        vol = permute(vol, [3,2,1]);
    case 'sagittal'
        vol = permute(vol, [1,3,2]);
    end

    %loop through the sub-cortical area index values and overlay them on the plot
    for ii=1:length(subCorticalInd)
        these_areas = subCorticalInd{ii};
        areaMask = any(vol==these_areas(1),3);

        if length(these_areas)>1
            for kk=2:length(these_areas)
                areaMask = areaMask+any(vol==these_areas(kk),3);
            end
        end
        areaMask = any(areaMask,3);

        hold on
        B=bwboundaries(areaMask);
        for kk=1:length(B)
            S(end+1)=plot(ax,B{kk}(:,2),B{kk}(:,1),'b--');
        end
        hold off
    end



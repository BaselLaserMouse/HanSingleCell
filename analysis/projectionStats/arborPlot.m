function varargout = arborPlot(laminarData,cleanCells,cleanOnly,indToPlot,stats)
    % Plot arbors in a target area
    %
    % out = arborPlot(laminarData,cleanCells,cleanOnly)
    %
    %
    % Inputs
    % laminarData - see below
    % cleanCells - xylem structure
    % cleanOnly - false by default



    if nargin<3 || isempty(cleanOnly)
        cleanOnly=false;
    end

    if cleanOnly %TODO: does not work -- causes bugs down-stream
        % Filter neurons by premature/not 
        s=scc.diagnostic.summariseTerminationTypes(cleanCells,true);
        s=s(4); %No back-labelled cells V1 only
        cleanCells = s.IDsofCleanCells;

        %Go through everything and remove non-clean cells
        for ii=1:length(laminarData)
            tData=laminarData(ii);
            ind=[];
            for jj=1:length(tData.cellID)
                if ~isempty( strmatch(tData.cellID{jj},cleanCells) )
                    ind(end+1)=jj;
                end
            end % jj
            tData.cellID = tData.cellID(ind);
            tData.counts = tData.counts(:,:,ind);
            laminarData(ii) = tData;
        end
    end





    %doDotPlot



        clf
        tData=laminarData(indToPlot);

        tPlanes = tData.planesInARA;
        tCells = tData.cellID;
        tAreaID = name2structureID(tData.areaName);

        c=cleanCells.returnData;

        details = [c.details];
        cellIDs = {details.cellID};

        subplot(2,2,1)
        hold on
        for ii=1:length(tCells)
            ind=strmatch(tCells{ii},cellIDs);
            plotTree(c(ind),tData,stats,2);
        end

        title(tData.areaName)
        axis ij equal
        drawnow 


        subplot(2,2,2)
        hold on
        for ii=1:length(tCells)
            ind=strmatch(tCells{ii},cellIDs);
            plotTree(c(ind),tData,stats,1);
        end

        title(tData.areaName)
        axis ij equal
        drawnow


        subplot(2,2,3)
        hold on
        for ii=1:length(tCells)
            ind=strmatch(tCells{ii},cellIDs);
            plotTree(c(ind),tData,stats,-1);
        end

        view(3)
        axis ij equal
        grid on
        box on




%-------------------------------------------------------------------------------------------
    function h=plotTree(data,tData,stats,dimPlot)
        if nargin<3
            stats=[];
        end
        if nargin<4
            dimPlot=2;
        end

        neuriteTree = data.neuriteTree;
        segmentRows = neuriteTree.getsegments;

        thisArea = name2structureID(tData.areaName);
        tColor=[0,0,0.15,0.3];

        for ii=1:length(segmentRows)
            theseNodes = neuriteTree.Node(segmentRows{ii});
            theseData = ones(length(theseNodes),3);


            theseInds = data.pointsInARA.rawSparseData.ARAindex(segmentRows{ii});
            for jj=1:length(theseNodes)
                if theseInds(jj)==thisArea
                    theseData(jj,:) = theseNodes{jj};
                else
                    theseData(jj,:) = nan;
                end
            end

            %Only keep planes which have axon in LM 
            %f=find(~ismember(round(theseData(:,1)), tData.planesInARA));
            %theseData(f,:) = nan;

            if ~isempty(stats)
                theseData(:,2) = theseData(:,2) - stats.mu.ml;
                theseData(:,3) = theseData(:,3) - stats.mu.dv;
                theseData(:,1) = theseData(:,1) - stats.mu.rc;
                fittedValues = stats.fitfunc(theseData(:,2),theseData(:,1), stats.b);
                theseData(:,3) = theseData(:,3) - fittedValues;
            end

            if dimPlot>0
                h(ii)=plot(theseData(:,dimPlot), theseData(:,3), '-', ...
                    'linewidth',0.5, 'Color', tColor);
            else
                h(ii)=plot3(theseData(:,1), theseData(:,2), theseData(:,3), '-', ...
                    'linewidth',0.5, 'Color', tColor);
                if ~isempty(stats)
                    h(ii)=plot3(theseData(:,1), theseData(:,2), fittedValues, '.r', ...
                    'linewidth',0.5, 'Color', tColor);
                end
            end
        end

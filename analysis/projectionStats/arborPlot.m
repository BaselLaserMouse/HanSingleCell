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
    % indToPlot -  which index of laminarData to plot
    % stats - the areaFits (see fitAreaPlane)
    %
    % e.g.
    % >> ii=3;arborPlot(lamData,cleanCells,false,ii,areaFits{ii});




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
        colors = parula(length(tCells));

        subplot(2,2,1)
        hold on
        for ii=1:length(tCells)
            ind=strmatch(tCells{ii},cellIDs);
            plotTree(c(ind),tData,stats,2,colors(ii,:));
        end

        title([tData.areaName, ' (coronal)'])
        axis ij equal


        %try adding layer boundaries
        addLayerBoundaries(stats,2)
        drawnow


        subplot(2,2,2)
        hold on
        for ii=1:length(tCells)
            ind=strmatch(tCells{ii},cellIDs);
            plotTree(c(ind),tData,stats,1,colors(ii,:));
        end

        title([tData.areaName, ' (sagittal)'])
        axis ij equal
        drawnow
        %try adding layer boundaries
        addLayerBoundaries(stats,1)
return 
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

        %try adding layer boundaries
        addLayerBoundaries(stats,-1)



%-------------------------------------------------------------------------------------------
    function h=plotTree(data,tData,stats,dimPlot,tColor)
        if nargin<3
            stats=[];
        end
        if nargin<4
            dimPlot=2;
        end
        if nargin<5
            tColor=[0,0,0.15,0.3];
        end

        neuriteTree = data.neuriteTree;
        segmentRows = neuriteTree.getsegments;

        thisArea = name2structureID(tData.areaName);


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
                theseData = applyTranform(theseData,stats);
            end

            if dimPlot>0
                h(ii)=plot(theseData(:,dimPlot), theseData(:,3), '-', ...
                    'linewidth',0.5, 'Color', tColor);
            else
                h(ii)=plot3(theseData(:,1), theseData(:,2), theseData(:,3), '-', ...
                    'linewidth',0.5, 'Color', tColor);
                if ~isempty(stats)
                %    h(ii)=plot3(theseData(:,1), theseData(:,2), fittedValues, '.r', ...
                 %   'linewidth',0.5, 'Color', tColor);
                end
            end
        end % for ii=...

    function addLayerBoundaries(stats,dimPlot)
        if ~isfield(stats,'layers')
            return
        end

        p=parula(length(stats.layers));
        minX=zeros(1,length(stats.layers));
        muY=zeros(1,length(stats.layers));
        for ii=1:length(stats.layers)
            tData = [stats.layers(ii).RC,stats.layers(ii).ML,stats.layers(ii).DV];
            tData = applyTranform(tData,stats);

            % Fit a surface to these layer data
            X = [ones(size(tData,1),1), tData(:,1), tData(:,2)];
            fitfunc = @(x,y,b) b(1) + b(2)*x + b(3)*y;
            [b,fSTATS] = regress(tData(:,3),X);

            x1fit = 1.2*min(tData(:,1)):1:max(tData(:,1))*1.2;
            x2fit = 1.2*min(tData(:,2)):1:max(tData(:,2))*1.2;
            [X1FIT,X2FIT] = meshgrid(x1fit,x2fit);
            if isempty(X1FIT)
                fprintf('Failed to generate meshgrid for layer %d\n', ii);
                continue
            end
            YFIT = fitfunc(X1FIT, X2FIT, b);

            if dimPlot>0
              %  plot(tData(:,dimPlot), tData(:,3),'.','color',p(ii,:))
                if dimPlot==1
                    n=round(size(X1FIT,1)/2);
                    X = X1FIT(n,:);
                    Y = YFIT(n,:);
                elseif dimPlot==2
                    n=round(size(X2FIT,2)/2);
                    X = X2FIT(:,n);
                    Y = YFIT(:,n);
                end
                plot(X,Y,':','color',[1,0,0,0.5],'linewidth',1)
                minX(ii)=min(X);
                muY(ii)=mean(Y);
            else
                plot3(tData(:,1), tData(:,2), tData(:,3), '.','color',p(ii,:));
            end
        end

        if dimPlot>0 %Add a scale bar and layer labels

            labX=median(minX);
            labelText={'L1','L2/3','L4','L5','L6'};
            for ii=1:length(minX)
                if muY(ii) == 0 && minX(ii)==0, continue, end
                if ii<length(minX)
                    labY = mean(muY(ii:ii+1));
                    text(labX,labY,labelText{ii},'FontSize',12)
                    plot([labX,labX+4], [muY(ii),muY(ii)], '--k')
                end
                plot([labX,labX+4], [muY(ii),muY(ii)], '--k')
            end

            %Scale bar
            xl = xlim;
            rightPoint = xl(2)-range(xl)*0.2;
            yl = ylim;
            ypos = yl(2)-range(yl)*0.05;

            plot([rightPoint-8,rightPoint], [ypos,ypos], '-k', 'LineWidth', 5)
            axis off
        end

    function transformedValues = applyTranform(points,stats)
        %Subtract the offset
        points(:,2) = points(:,2) - stats.mu.ml;
        points(:,3) = points(:,3) - stats.mu.dv;
        points(:,1) = points(:,1) - stats.mu.rc;

        %Apply the affine tranform (ordered ml, rc, dv)
        tPoints = points(:,[2,1,3])';
        affineTransformedPoints = applyAffineTransform(tPoints,stats.affine,[],false);

        transformedValues = affineTransformedPoints([2,1,3],:)'; %Flip back the dimensions to keep things consistent here

        % Now subtract the curviture of the surface from the DV values
        fittedValues = stats.fitfunc(transformedValues(:,2), transformedValues(:,1), stats.b);
        transformedValues(:,3) = transformedValues(:,3)-fittedValues;



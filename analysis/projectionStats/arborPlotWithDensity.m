function varargout = arborPlotWithDensity(laminarData,cleanCells,cleanOnly,indToPlot,stats)
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

        ax1=axes('Position', [0.02,0.1,0.30,0.80]);
        hold on
        treeData=[];
        for ii=1:length(tCells)
            ind=strmatch(tCells{ii},cellIDs);
            [~,tmp]=plotTree(c(ind),tData,stats,2,colors(ii,:));
            treeData=[treeData;tmp];
        end

        title([tData.areaName, ' (coronal)'])

        set(ax1,'XLimMode','manual');

        addLayerBoundaries(stats,2)
        YL=ylim;
        axis equal off

        p=plotboxpos(ax1);
        axes('Position', [0.32,p(2),0.10,p(4)])
        plotDensity(treeData)
        set(gca,'XLim',YL)
        drawnow
        box on


        %Sagittal
        ax3=axes('Position', [0.52,0.1,0.30,0.80]);
        hold on
        for ii=1:length(tCells)
            ind=strmatch(tCells{ii},cellIDs);
            plotTree(c(ind),tData,stats,1,colors(ii,:));
        end

        title([tData.areaName, ' (sagittal)'])
        set(ax3,'XLimMode','manual')
        
        %try adding layer boundaries
        addLayerBoundaries(stats,1)
        axis equal off
        p=plotboxpos(ax3);
        axes('Position', [0.82,p(2),0.10,p(4)])
        plotDensity(treeData)
        set(gca,'XLim',YL)




%-------------------------------------------------------------------------------------------
    function [h,treeData]=plotTree(data,tData,stats,dimPlot,tColor)
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

        treeData=[];
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

            h(ii)=plot(theseData(:,dimPlot), theseData(:,3), '-', ...
                'linewidth',0.5, 'Color', tColor);
            treeData=[treeData;theseData(:,3)];
        end % for ii=...
        treeData(isnan(treeData))=[];
        axis ij equal tight

    function plotDensity(treeData)
        [n,x]=hist(treeData,0:1:round(max(treeData)));
        n = smooth(n,3);
        n(1)=0; %Otherwise the plot baseline is skewed

        ptch=patch(x,n,1);

        set(gca,'view',[90 90])
        set(ptch,'FaceColor','k','EdgeColor','k')
        axis off

    function addLayerBoundaries(stats,dimPlot)
        if ~isfield(stats,'layers')
            return
        end

        layerTickSize=3; %Size of the layer ticks in 25 micron voxels;
        minX=zeros(1,length(stats.layers));
        maxX=zeros(1,length(stats.layers));
        muY=zeros(1,length(stats.layers));
        for ii=1:length(stats.layers)
            tData = [stats.layers(ii).RC,stats.layers(ii).ML,stats.layers(ii).DV];
            tData = applyTranform(tData,stats);

            % Fit a surface to these layer data
            X = [ones(size(tData,1),1), tData(:,1), tData(:,2)];
            fitfunc = @(x,y,b) b(1) + b(2)*x + b(3)*y;
            [b,fSTATS] = regress(tData(:,3),X);

            x1fit = min(tData(:,1))-layerTickSize:1:max(tData(:,1))+layerTickSize;
            x2fit = min(tData(:,2))-layerTickSize:1:max(tData(:,2))+layerTickSize;
            [X1FIT,X2FIT] = meshgrid(x1fit,x2fit);
            if isempty(X1FIT)
                fprintf('Failed to generate meshgrid for layer %d\n', ii);
                continue
            end
            YFIT = fitfunc(X1FIT, X2FIT, b);

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
            plot(X,Y,':','color',[1,0,0,0.7],'linewidth',1)
            minX(ii)=min(X);
            maxX(ii)=max(X);
            muY(ii)=mean(Y);

        end %for ii=1:length(stats.layers)

        %The little things on the edges
        labX=xlim;

        labelText={'L1','L2/3','L4','L5','L6'};
        for ii=1:length(minX)
            if muY(ii) == 0 && minX(ii)==0, continue, end
                if ii<length(minX)
                    labY = mean(muY(ii:ii+1));
                    text(labX(1),labY,labelText{ii},'FontSize',12)
                end
            plot([labX(1),labX(1)+layerTickSize], [muY(ii),muY(ii)], '--k')
            plot([labX(2),labX(2)-layerTickSize], [muY(ii),muY(ii)], '--k')
        end

        %Scale bar
        xl = xlim;
        rightPoint = xl(2)-range(xl)*0.2;
        yl = ylim;
        ypos = muY(end)+3;

        plot([rightPoint-8,rightPoint], [ypos,ypos], '-k', 'LineWidth', 5)
        xlim(labX)

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



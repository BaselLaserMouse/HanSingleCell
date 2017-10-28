function varargout = arborPlot_3D(laminarData,cleanCells,cleanOnly,indToPlot,stats)
    % Plot arbors in a target area
    %
    % out = arborPlot_3D(laminarData,cleanCells,cleanOnly)
    %
    % Purpose
    % Used for brain areas with a lot curviture 
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
    % >> ii=3;arborPlot-3D(lamData,cleanCells,false,ii,areaFits{ii});




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




    clf
    tData=laminarData(indToPlot);

    tPlanes = tData.planesInARA;
    tCells = tData.cellID;
    tAreaID = name2structureID(tData.areaName);

    c=cleanCells.returnData;

    details = [c.details];
    cellIDs = {details.cellID};
    colors = parula(length(tCells));

    hold on
    s=stats;
    %s.fDims = 1:3;
    for ii=1:length(tCells)
        ind=strmatch(tCells{ii},cellIDs);
        plotTree(c(ind),tData,s);
    end

    layerPlotter(stats)

    view(3)
    axis ij equal
    grid on
    box on




%-------------------------------------------------------------------------------------------
    function h=plotTree(data,tData,stats)
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

            h(ii)=plot3(theseData(:,1), theseData(:,2), theseData(:,3), '-', ...
                'linewidth',0.5, 'Color', tColor);

        end % for ii=...

    function layerPlotter(stats)
        % Plots layers and fits them
        p=parula(length(stats.layers));



        for ii=1:length(stats.layers)
            tData = [stats.layers(ii).RC,stats.layers(ii).ML,stats.layers(ii).DV];
            tData = applyTranform(tData,stats);


            % Fit a surface to these data
            X = [ones(size(tData,1),1), tData(:,1), tData(:,2)];
            fitfunc = @(x,y,b) b(1) + b(2)*x + b(3)*y;
            [b,fSTATS] = regress(tData(:,3),X);

            x1fit = 1.2*min(tData(:,1)):1:max(tData(:,1))*1.2;
            x2fit = 1.2*min(tData(:,2)):1:max(tData(:,2))*1.2;
            [X1FIT,X2FIT] = meshgrid(x1fit,x2fit);
            YFIT = fitfunc(X1FIT, X2FIT, b);

            %Fit it
            %X = [ones(size(tData,1),1), tData(:,1), tData(:,2), tData(:,1).^2, tData(:,2).^2 tData(:,1).^3, tData(:,2).^3];
            %X = [ones(size(tData,1),1), tData(:,1), tData(:,2), tData(:,1).^2, tData(:,2).^2];


            %fitfunc = @(x,y,b) b(1) + b(2)*x + b(3)*y + b(4)*x.^2 + b(5)*y.^2 + b(6)*x.^3 + b(7)*y.^3;
            %fitfunc = @(x,y,b) b(1) + b(2)*x + b(3)*y + b(4)*x.^2 + b(5)*y.^2;

            plot3(tData(:,1), tData(:,2), tData(:,3), '.','color',p(ii,:));

        end %for




    function transformedValues = applyTranform(points,stats)
        %Subtract the offset 

        % This way of handling the dimensions is crappy, but it works. 
        if all(stats.fDims==1:3) %No dim flipping
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

        elseif all(stats.fDims==[3,2,1])
            % This flip is for a swap of RC and DV
            points(:,2) = points(:,2) - stats.mu.ml;
            points(:,3) = points(:,3) - stats.mu.rc; %replace .dv with .rc 
            points(:,1) = points(:,1) - stats.mu.dv; %replace .rc with .dv

            %Apply the affine tranform (should be ordered ml, dv, rc)
            tPoints = points(:,[2,3,1])';
            affineTransformedPoints = applyAffineTransform(tPoints,stats.affine,[],false);

            transformedValues = affineTransformedPoints([2,3,1],:)'; %Flip back the dimensions to keep things consistent here            

        elseif all(stats.fDims==[2,1,3])
            % This flip is for a swap of ML and RC
            points(:,2) = points(:,2) - stats.mu.rc; %replace ml with rc
            points(:,3) = points(:,3) - stats.mu.dv;
            points(:,1) = points(:,1) - stats.mu.ml; %replace rc with ml

            %Apply the affine tranform (ordered ml, rc, dv)
            tPoints = points(:,[1,2,3])';
            affineTransformedPoints = applyAffineTransform(tPoints,stats.affine,[],false);

            transformedValues = affineTransformedPoints([1,2,3],:)'; %Flip back the dimensions to keep things consistent here   

        else
            error('Unknown dim flip')
        end





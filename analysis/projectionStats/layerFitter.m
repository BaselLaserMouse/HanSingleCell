function varargout = layerFitter(stats)
    % function varargout = layerFitter(stats)
    %
    % viewer for fitAreaPlane




    if ~isfield(stats,'layers')
        return
    end
    clf

    subplot(2,2,1)
    hold on
    layerPlotter(stats,2)
    title(stats.areaName)

    subplot(2,2,2)
    hold on
    layerPlotter(stats,1)

    subplot(2,2,3)
    hold on
    layerPlotter(stats,-1)
    view(3)




function layerPlotter(stats,dimPlot)
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

        if dimPlot>0
            plot(tData(:,dimPlot), tData(:,3),'.','color',p(ii,:))
            if dimPlot==1
                n=round(size(X1FIT,1)/2);
                X = X1FIT(n,:);
                Y = YFIT(n,:);
            elseif dimPlot==2
                n=round(size(X2FIT,1)/2);
                X = X2FIT(:,n);
                Y = YFIT(:,n);
            end
            plot(X,Y,'--','color',[1,1,1]*0.25,'linewidth',2)

        else
            plot3(tData(:,1), tData(:,2), tData(:,3), '.','color',p(ii,:));

            M=mesh(X1FIT,X2FIT,YFIT);
            M.FaceAlpha=0.5;

        end
    end


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

            %NO CURVE CORRECTION
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



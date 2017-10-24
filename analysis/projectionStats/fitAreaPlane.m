function varargout = fitAreaPlane(areaInd,atlasVol,E)
    % Fit a plane to a named brain area
    %
    % function fitStats = fitAreaPlane(areaInd,atlasVol)
    %
    % Purpose
    % We want to make nice laminar plots of axons in a particular brain area, however
    % the brain surface is curved. So we want to be able to "look" in a direction
    % orthogonal to the cortical layers. We also want to subtract any local curvature. 
    % This function performs these operations using the following steps:
    % 1) Find the geometric centre of L1 of the cortical area in question.
    % 2) Make the area parallel to the plot axes by rotating it around this point using 
    %    an affine transformation. 
    % 3) Fit a 2-D curve to the surface to remove gross local distortion. 
    % 
    % These operations can be applied to points, but they will need to be applied in order.
    % 
    %
    % Inputs
    % areaInd - The brain area to fit. Can be one of the following:
    %          a) Name of a cortical area. e.g. `Primary visual area'
    %          b) Index of a cortical area. e.g. 385  (which is V1)
    %          c) Index of a layer in a cortical area. e.g. 33 (which is V1 L6a)
    %
    %          If the user supplies (c), the surface of this is fitted and the parameters
    %          returned. If the user supplies (a) or (b) then the functions finds L1 and 
    %          fits this.... TODO: fit other layer too (as offsets?).
    % atlasVol - optional. By default the ARA is pulled in using aratools.atlascacher.getCachedAtlas
    %          You can also supply either this structure or the array manually using this argument.
    %          This should be a bit faster.
    % E - canny edge detected version of the atlas. do: E=canny(atlasVol);
    %
    %
    %  TODO: decide how to handle the hemispheres: fit both or just one? They should be the same
    %        with opposite signs for the ML slope parameters...


    if nargin<1 || isempty(areaInd)
        areaInd = 421;
    else
        %otherwise handle the input argument
        if ischar(areaInd)
            areaInd = name2structureID(areaInd);
        end
        childTable = getAllenStructureList('childrenOf',areaInd);

        if size(childTable,1)>1
            %Find L1 or quit
            f=find(...
                cell2mat(...
                    cellfun(@(x) ~isempty(regexp(x,'.*, [Ll]ayer 1')) , childTable.name, 'UniformOutput', false)...
                    )...
                );
            if isempty(f)
                fprintf('Can not find "layer 1" in area %d\n', areaInd);
                return
            elseif length(f)>1
                fprintf('Found multiple "layer 1" in area %d\n', areaInd);
                return
            end

            fprintf('Fitting %s\n', childTable.name{f});
            areaInd = childTable.id(f);


        end % close size(childTable,1)>1

    end % nargin<1 || isempty(areaInd)


    if nargin<2
        a=aratools.atlascacher.getCachedAtlas;
        atlasVol = a.atlasVolume;
    end

    if isstruct(atlasVol)
        atlasVol = atlasVol.atlasVolume;
    end



    %Work only on the hemisphere where we have acquired data
    atlasVol(:,1:round(size(atlasVol,2)/2),:)=0;
    E(:,1:round(size(atlasVol,2)/2),:)=0;

    % DIM FLIP
    %atlasVol = permute(atlasVol,[3,2,1]);
    %E = permute(E,[3,2,1]);

    %Find the surface points
    ef=find(E==1);
    [DVe,MLe,RCe]=ind2sub(size(E),ef(1:5:end));



    mask = (atlasVol==areaInd); 


    %Now get the coordinates of these points
    f=find(mask+E == 2);
    %DV: dorso-ventral, ML: medio-lateral, RC: rostrocaudal
    [DV,ML,RC] = getAreaSurfacePixels(atlasVol,E,areaInd);

    %Get the coords of all the layers
    for ii=1:length(childTable.name)
        layerCoords(ii).name = childTable.name{ii};
        layerCoords(ii).id = childTable.id(ii);
        [layerCoords(ii).DV, layerCoords(ii).ML, layerCoords(ii).RC] = getAreaSurfacePixels(atlasVol,E,childTable.id(ii));
    end


    clf
    subplot(2,2,1)
    plot3(MLe,RCe,DVe,'.k') %whole brain

    hold on
    plot3(ML,RC,DV,'.', 'Color', [1,0,0.33]) %Selected area

    xlabel('ML'), ylabel('RC'), zlabel('DV')
    axis equal 
    hold on
    plot3(mean(ML),  mean(RC), mean(DV), 'or', 'MarkerFaceColor', 'r')
    axis ij, box on, grid on
    title(['Estimated surface of area ', structureID2name(areaInd)])
    %ylim([100,300])

    %Add a plane to indicate the surface we try to register to
    xl=xlim;
    yl=ylim;
    ptch=patch([xl(1),xl(2),xl(2),xl(1)], ....
        [yl(1),yl(1),yl(2),yl(2)], 1);
    set(ptch,'FaceColor','b','FaceAlpha',0.25)

    %Now we centre around the mean and we'll perform all transformations around this
    dv = DV-mean(DV);
    ml = ML-mean(ML);
    rc = RC-mean(RC);

    dvE = DVe-mean(DV);
    mlE = MLe-mean(ML);
    rcE = RCe-mean(RC);


    % DV is the dorso-ventral position that we want to "straighten"
    % Now we make the appropriate matrices to feed absor
    orig = [ml,rc,dv]';
    origE = [mlE,rcE,dvE]';
    target = orig; target(3,:)=0;
    affineStats = absor(orig,target);


    % Apply the transformation
    subplot(2,2,2)
    affineTransformed = applyAffineTransform(orig,affineStats,target);
    affineTransformedE = applyAffineTransform(origE,affineStats,[],false);
    title('Result of affine transform')

    % Now we attempt to correct for surface curvature 
    mlAT = affineTransformed(1,:)';
    rcAT = affineTransformed(2,:)';
    dvAT = affineTransformed(3,:)';

    mlATE = affineTransformedE(1,:)';
    rcATE = affineTransformedE(2,:)';
    dvATE = affineTransformedE(3,:)';

    X = [ones(size(DV)) mlAT rcAT mlAT.^2 rcAT.^2 mlAT.^3 rcAT.^3];
    X = [ones(size(DV)) mlAT rcAT mlAT.^2 rcAT.^2];
    [b,BINT,R] = regress(dvAT,X);



    fitfunc = @(x,y,b) b(1) + b(2)*x + b(3)*y + b(4)*x.^2 + b(5)*y.^2 + b(6)*x.^3 + b(7)*y.^3;
    fitfunc = @(x,y,b) b(1) + b(2)*x + b(3)*y + b(4)*x.^2 + b(5)*y.^2;

    % The fitted surface
    subplot(2,2,3)
    plot3(mlATE,rcATE,dvATE,'.k')
    hold on
    plot3(mlAT,rcAT,dvAT,'.r')


    x1fit = 1.5*min(ml):1:max(ml)*1.5;
    x2fit = 1.5*min(rc):1:max(rc)*1.5;
    [X1FIT,X2FIT] = meshgrid(x1fit,x2fit);

    YFIT = fitfunc(X1FIT, X2FIT, b);
    M=mesh(X1FIT,X2FIT,YFIT);
    xlabel('ML'), ylabel('RC'), zlabel('DV')
    axis equal
    title('Fit to brain area surface')
    %ylim([-75,100])
    zlim([-20,100])

    % Subtract the fit from the surface (the residuals if we have just one layer)
    subplot(2,2,4)

    DVfit = fitfunc(mlAT,rcAT,b);

    plot3(mlAT,rcAT,dvAT-DVfit,'r.')
    xlabel('ML'), ylabel('RC'), zlabel('DV')
    axis equal
    box on, grid on

    title('Fit to brain area surface')


    if nargout>0
        stats.b = b;
        stats.r = R;
        stats.mu.dv = mean(DV);
        stats.mu.ml = mean(ML);
        stats.mu.rc = mean(RC);
        stats.fitfunc = fitfunc;
        stats.affine = affineStats;
        stats.areaInd = areaInd;
        stats.areaName = structureID2name(areaInd);
        stats.layers = layerCoords;
        varargout{1}=stats;

    end




function [DV,ML,RC] = getAreaSurfacePixels(atlasVol,E,areaInd)
    %DV: dorso-ventral
    %ML: medio-lateral
    %RC: rostrocaudal

    mask = (atlasVol==areaInd); 
    %Now get the coordinates of these points
    f=find(mask+E == 2);
    [DV,ML,RC]=ind2sub(size(mask),f);
    %DV: dorso-ventral
    %ML: medio-lateral
    %RC: rostrocaudal


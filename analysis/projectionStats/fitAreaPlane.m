function varargout = fitAreaPlane(areaInd,atlasVol)
    % Fit a plane to a named brain area
    %
    % function fitStats = fitAreaPlane(areaInd,atlasVol)
    %
    % Purpose
    % Fit plane to brain area surface and return fit paramaters.
    % This can be used to "straighten" an area for plotting.
    % 
    %
    % Inputs
    % areaInd - The brain area to fit. Can be one of the following:
    %          a) Name of a cortical area. e.g. `Primary visual area'
    %          b) Index of a cortical area. e.g. 385  (which is V1)
    %          c) Index of a layer in a cortica area. e.g. 33 (which is V1 L6a)
    %
    %          If the user supplies (c), the surface of this is fitted and the parameters
    %          returned. If the user supplies (a) or (b) then the functions finds L1 and 
    %          fits this.... TODO: fit other layer too (as offsets?).
    % atlasVol - optional. By default the ARA is pulled in using aratools.atlascacher.getCachedAtlas
    %          You can also supply either this structure or the array manually using this argument.
    %          This should be a bit faster.
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



    %Work on one hemisphere only 
    atlasVol(:,1:round(size(atlasVol,2)/2),:)=[];

    mask = (atlasVol==areaInd);



    % Keep the surface of the area
    d=mask-circshift(mask,1,1);
    d(d<1)=0;

    %Now get the coordinates of these points
    f=find(d>0);
    [DV,ML,RC]=ind2sub(size(d),f);

    %DV - dorso-ventral
    %ML - medio-lateral
    %RC - rostrocaudal


    clf

    %subplot(2,2,1), imagesc(squeeze( sum(d,1)) )
    %subplot(2,2,2), imagesc(squeeze( sum(d,3)) )
    

    % The extracted surface points
    subplot(2,2,1)

    plot3(ML,RC,DV,'.', 'Color', [0,0,0.33])
    xlabel('ML'), ylabel('RC'), zlabel('DV')
    axis equal 
    hold on
    plot3(mean(ML),  mean(RC), mean(DV), 'or', 'MarkerFaceColor', 'r')
    axis ij, box on, grid on
    title(['Estimated surface of area ', structureID2name(areaInd)])

    %Now we fit a surface around the mean by predicting the DV position based on ML and RC
    %dv = DV-mean(DV);
    %ml = ML-mean(ML);
    %rc = RC-mean(RC);


    X = [ones(size(DV)) ML RC ML.^2 RC.^2 ML.^3 RC.^3];
    [b,BINT,R] = regress(DV,X);


    fitfunc = @(x,y,b) b(1) + b(2)*x + b(3)*y + b(4)*x.^2 + b(5)*y.^2 + b(6)*x.^3 + b(7)*y.^3;

    % The fitted surface
    subplot(2,2,2)
    scatter3(ML,RC,DV,'filled')
    hold on
    x1fit = min(ML):1:max(ML);
    x2fit = min(RC):1:max(RC);
    [X1FIT,X2FIT] = meshgrid(x1fit,x2fit);
    YFIT = b(1) + b(2)*X1FIT + b(3)*X2FIT + b(4)*X1FIT.^2 + b(5)*X2FIT.^2 + b(6)*X1FIT.^3 + b(7)*X2FIT.^3;
    YFIT = fitfunc(X1FIT, X2FIT, b);
    M=mesh(X1FIT,X2FIT,YFIT);
    xlabel('ML'), ylabel('RC'), zlabel('DV')
    axis equal
    title('Fit to brain area surface')

    % Subtract the fit from the surface (the residuals if we have just one layer)
    subplot(2,2,3)

    DVfit = fitfunc(ML,RC,b);% b(1) + b(2)*ML + b(3)*RC + b(4)*ML.^2 + b(5)*RC.^2 + b(6)*ML.^3 + b(7)*RC.^3;

    plot3(ML,RC,DV-DVfit,'r*')
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
        varargout{1}=stats;
    end

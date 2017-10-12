function varargout = laminarPlots(cellMat,cleanCells,areaMap)
    % Plot distribution of terminals by layer within a subset of defined areas
    %
    % varargout = laminarPlots(cellMat,cleanCells,areaMap)
    %
    %
    % Inputs
    % cellMat - the binary matrices and associated labels from Justus packages in structure.
    %           if the structure is of length 2 then we assume the first index contains
    %           the "good" cells from the main paper and the second index the premature
    %           cells.
    % cleanCells - Xylem data structure containing the selected good cells (e.g. no crappy 
    %              ones where the brain was falling apart)
    % areaMap - is made by makeAreaMap
    %
    % >> load ~/tvtoucan/Mrsic-Flogel/hanyu/Analyses/cleanCells.mat
    % >> load allCellMat.mat
    % >> load areaMap
    % >> projectionDensity(allCellMat,cleanCells,m)


    % Get the length of axon for each neuron in the cellMat array

    cellMat(1).totalDistance=[];
    for ii=1:length(cellMat)
        cellMat(ii) = getAxonLength(cellMat(ii), cleanCells);
    end


    clf
    [~,t]=aratools.utils.whiteMatterInds;
    D=pointsByAreaPlot(cleanCells,'dataMetric', 'upSampledPoints', ...
        'excludeBorders', 4, ...
        'excludeAreas', {'Intercalated amygdalar nucleus','Out of brain','lateral ventricle','Primary visual','Cortical subplate','ventricular',t{:}});

    M=D.dataMat*5; % it's 5 microns per voxel
    M(M<20E2)=0; % Exclude areas with small lengths
    M = M * 1E-3; % convert to mm
    % TODO: also check what "total distance" is in the cleanCells array. What is the distance between points?

    %find the columns in the matrix that correspond to the premature cells so we can highight these differently
    prem = zeros(1,size(M,2),'logical');
    for ii=1:length(prem)
        if ~isempty(strmatch(D.cellIDs{ii},cellMat(2).id))
            prem(ii) = true;
        end
    end


    % Axon length as a function of the number of target areas
    subplot(2,2,1)
    p = M(:,prem);
    plot(sum(p>0), sum(p),'sb')
    hold on
    n = M(:,~prem);
    plot(sum(n>0), sum(n),'sb');%ob','markerfacecolor',[0.5,0.5,1])
    hold off
    ylabel('Axon length in all target areas (mm)')
    xlabel('Number of projection targets')
    xlim([0.5 7])
    jitter
    title('Total axon length across all target areas as function of the number of target areas')
    grid 
    box off



    % To calculate axon density averaged over all areas as a function of number of target areas
    % we will first make a version of M where we don't have axon length in mm but axon length in mm 
    % divided by area volume
    Md = M;
    vols = zeros(length(D.areaNamesInSamples),1);
    for ii=1:length(vols)
        tArea = D.areaNamesInSamples{ii};
        if ~areaMap.isKey(tArea)
            fprintf('No key %s found in map. skipping\n', tArea)
            continue
        end
        vols(ii) = areaMap(tArea);
        if vols(ii)==0;
            fprintf('Warning: returning 0 mm^3 for area %s\n', tArea)
        end
    end
    vols = repmat(vols,1,size(M,2));
    Md = M./vols;
    Md(Md==0) = NaN ;
    Md(isinf(Md)) = NaN ;



    subplot(2,2,2)
    mu = nanmean(Md,1);
    s = sum(~isnan(Md),1);
    plot(s,mu,'ob');
    xlim([0.5 7])
    jitter 
    xlabel('Number of projection targets')
    ylabel('Average axon length/area size in mm')
    title('Axon density averaged over all areas as function of num target areas')
    grid 
    box off



    % length as a function of axon density
    subplot(2,2,3)
    cla 
    hold on
    plot(Md(:),M(:),'.','Color',[1,1,1]*0.5)
    for ii=1:size(Md,1)
        tMd = Md(ii,:);
        tM = M(ii,:);

        for kk=1:length(tM)
            text(tMd(kk), tM(kk), num2str(ii))
        end

    end

    title(sprintf('9. %s; 3. %s; 28. %s\n11. %s 12. %s', ...
        D.areaNamesInSamples{9}, ...
        D.areaNamesInSamples{3}, ...
        D.areaNamesInSamples{28}, ...
        D.areaNamesInSamples{11}, ...
        D.areaNamesInSamples{11} ...
        ))


    grid 
    box off
    xlabel('Axon density (length/area)')
    ylabel('Axon length')



    % Does projection density to a chosen area vary depending on the number of target areas?
    subplot(2,2,4)
    areasToPlot = {'Lateral visual area', ...
        'posteromedial visual area', ...
        'Postrhinal area'};
    cla
    hold on
    mrkr = {{'or', 'MarkerFaceColor', [1.0,0.5,0.5]}, ...
            {'sb', 'MarkerFaceColor', [0.5,0.5,1.0]}, ...
            {'gp', 'MarkerFaceColor', [0.5,1.0,0.5]}};
    for ii=1:length(areasToPlot)
        ind=strmatch(areasToPlot{ii},D.areaNamesInSamples); % Row with connections to this are
        tmp = M(ind,:);
        f = find(tmp>0); % All neurons projecting to this area
        tmp = M(:,f); % The matrix with only these cells
        plot(sum(tmp>0), tmp(ind,:), mrkr{ii}{:});%ob','markerfacecolor',[0.5,0.5,1])
    end

    grid on
    xlim([0.5,5])
    jitter 

    legend(areasToPlot)
    ylabel('Axon length in named area (mm)')
    xlabel('Number of projection targets')
    title('Axon length in named target area as a function of the number of areas')

    if nargout>0
        varargout{1} = cellMat;
    end

function C = getAxonLength(C,cleanCells)
     D = cleanCells.returnData;
     details = [D.details];
     cellIDs = {details.cellID};

     for ii=1:length(C.id)
         ind = strmatch(C.id{ii}, cellIDs);
         C.totalDistance(ii) = D(ind).totalDistance;
     end
function varargout = projectionDensity(cleanCells,areaMap,plotPrem,plotNonV1)
    % Address whether the total length of axon scales with the number of target brain ares
    %
    % function varargout = projectionDensity(cellMat,cleanCells,areaMap)
    %
    %
    % Inputs
    % cleanCells - Xylem data structure containing the selected good cells (e.g. no crappy 
    %              ones where the brain was falling apart)
    % areaMap - is made by makeAreaMap
    % plotPrem - include premature termination cells
    % plotNonV1 - include cells with somata outside of V1
    %
    % e.g.
    % >> load ~/tvtoucan/Mrsic-Flogel/hanyu/Analyses/cleanCells.mat
    % >> load areaMap
    % >> projectionDensity(cleanCells,m)
    %

    if nargin<3
        plotPrem=false;
    end

    if nargin<4
        plotNonV1=false;
    end

    clf
    [~,t]=aratools.utils.whiteMatterInds;
    D=pointsByAreaPlot(cleanCells,'dataMetric', 'upSampledPoints', 'excludeBorders', 2, ...
        'excludeSomataNonV1', false, ...
        'excludeAreas', {'Intercalated amygdalar nucleus','Out of brain','lateral ventricle','Primary visual','Cortical subplate','ventricular systems',t{:}});
    M = D.dataMat;
    M(M<1)=0; % Exclude areas with small lengths

    S=scc.diagnostic.summariseTerminationTypes(cleanCells,true);

    Rdata = cleanCells.returnData;
    details = [Rdata.details];
    allCellIDs = {details.cellID};

    % Find the columns in the matrix that correspond to the premature cells so we can highlight these differently
    premIDs = allCellIDs(~S(4).cleanCell & S(4).allCells);
    premInd = zeros(1,length(D.cellIDs), 'logical');
    for ii=1:length(D.cellIDs)
        if ~isempty(strmatch(D.cellIDs{ii}, premIDs))
            premInd(ii) = true;
        end
    end

    %These are the clean V1 non-backlabelled neurons
    cleanV1ind = S(4).indexValuesOfcleanCells;
    cleanV1s = S(4).IDsofCleanCells;

    % Find the columns in the matrix that correspond to non-V1 cells so we can highlight these differently 
    nonV1ind = S(5).indexValuesOfAllCells;
    nonV1CellIDs = S(5).IDsOfAllCells;

    % Define marker styles
    cleanMrkr = {'ob', 'MarkerFaceColor', [0.5,0.5,1.0],'markersize',8};
    prematureMrkr = {'sb', 'MarkerFaceColor', [0.75,0.75,1.0],'markersize',8};
    nonV1mrkr = {'MarkerEdgeColor','r','markersize',9,'linewidth',2};

    % Axon length as a function of the number of target areas
    subplot(2,2,1)
    %Plot only clean v1 cells
    V1_plt_tmp = M(:,cleanV1ind);
    yV1=sum(V1_plt_tmp)
    xV1=sum(V1_plt_tmp>0)
    yV1./xV1
    nanmean(yV1./xV1)
    nanstd(yV1./xV1)
    plot(xV1, yV1, cleanMrkr{:});
    hold on 

    if plotPrem
        %Add premature cells
        tmp = M(:,premInd);
        plot(sum(tmp>0), sum(tmp), prematureMrkr{:});
    end


    %Add non-V1 cells (non-premature)
    if plotNonV1
        tmp = M(:,S(5).indexValuesOfcleanCells);
        fprintf('%d non-v1 premature non-premature cells\n', length(S(5).indexValuesOfcleanCells))
        plot(sum(tmp>0), sum(tmp), cleanMrkr{:}, nonV1mrkr{:});
    end


    %Add non-V1 cells (premature)
    if plotNonV1 && plotPrem
        tInd = setxor(S(5).indexValuesOfcleanCells,S(5).indexValuesOfAllCells);
        tmp = M(:,tInd);
        fprintf('%d non-v1 non-premature cells\n', length(tInd))
        plot(sum(tmp>0), sum(tmp), prematureMrkr{:}, nonV1mrkr{:});
    end

    jitter
    xlim([0.5,8.5])

    hold off
    ylabel('Axon length in all target areas (mm)')
    xlabel('Number of projection targets')
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
            fprintf('projectionDensity - No key %s found in map. skipping\n', tArea)
            continue
        end
        vols(ii) = areaMap(tArea);
        if vols(ii)==0;
            fprintf('Warning: returning a volume of 0 mm^3 for area "%s"\n', tArea)
        end
    end
    vols = repmat(vols,1,size(M,2));
    Md = M./vols;
    Md(Md==0) = NaN ;
    Md(isinf(Md)) = NaN ;



    % length as a function of axon density
    subplot(2,2,2)
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

    title(sprintf('16. %s; 13. %s; 5. %s\n18. %s 14. %s 15. %s', ...
        D.areaNamesInSamples{16}, ...
        D.areaNamesInSamples{13}, ...
        D.areaNamesInSamples{5}, ...
        D.areaNamesInSamples{18}, ...
        D.areaNamesInSamples{14}, ...
        D.areaNamesInSamples{15} ...
        ))


    grid 
    box off
    xlabel('Axon density (length/area)')
    ylabel('Axon length')



    subplot(2,2,3)
    mu = nanmean(Md,1);
    s = sum(~isnan(Md),1);
    plot(s,mu,'ok','MarkerFaceColor',[1,1,1]*0.5);
    xlim([0.5 9])
    jitter 
    xlabel('Number of projection targets')
    ylabel('Average axon length/area size in mm')
    title('Axon density averaged over all areas as function of num target areas (all cells)')
    grid 
    box off


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

    %Work with V1 only no backlabelled cells. Both premature and non-premature: 
    if plotPrem
        mTMP = M(:,S(4).indexValuesOfAllCells);
    else
        mTMP = M(:,S(4).indexValuesOfcleanCells);
    end
    for ii=1:length(areasToPlot)
        ind=strmatch(areasToPlot{ii},D.areaNamesInSamples); % Row with connections to this are
        tmp = mTMP(ind,:); % Connections to this area from V1 cells

        f = find(tmp>0); % All neurons projecting to this area
        tmp = mTMP(:,f); % The matrix with only these cells
        plot(sum(tmp>0), tmp(ind,:), mrkr{ii}{:});%ob','markerfacecolor',[0.5,0.5,1])
    end

    grid on
    xlim([0.5,6])
    jitter 

    legend(areasToPlot)
    ylabel('Axon length in named area (mm)')
    xlabel('Number of projection targets')
    title('Axon length in named target area as a function of the number of areas')

    if nargout>0
        varargout{1} = cellMat;
    end

    set(gcf,'PaperPosition',[0,0,40,40])
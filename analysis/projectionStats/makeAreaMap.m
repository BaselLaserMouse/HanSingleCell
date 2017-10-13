function m = makeAreaMap(cellMat,cleanCells)
    % Create a map structure of area volumes based on the unique names in allCellMat
    %
    % function m = makeAreaMap(cellMat,cleanCells)
    %
    %
    % Inputs
    % cellMat - the binary matrices and associated labels from Justus packages in structure.
    %           if the structure is of length 2 then we assume the first index contains
    %           the "good" cells from the main paper and the second index the premature
    %           cells.
    % cleanCells - the clean cells cell structure
    %
    %
    % Output
    % m - a map structure with area names as keys and mm^3 of area volume as values
    %
    % e.g.
    % >> load ~/tvtoucan/Mrsic-Flogel/hanyu/Analyses/cleanCells.mat
    % >> load allCellMat.mat
    % >> m = makeAreMap(allCellMat,cleanCells);
    %
    % Rob Campbell

    m = containers.Map;

    areaNames = unique(cat(1,cellMat.lab));
    [~,t]=aratools.utils.whiteMatterInds;
    D=pointsByAreaPlot(cleanCells,'dataMetric', 'upSampledPoints', 'excludeBorders', 2, ...
        'excludeSomataNonV1', false, ...
        'excludeAreas', {'Intercalated amygdalar nucleus','Out of brain','lateral ventricle','Primary visual','Cortical subplate','ventricular systems',t{:}});

    areaNames = unique([areaNames; D.areaNamesInSamples]);



    for ii=1:length(areaNames)
        vol = aratools.utils.getAreaVolume(areaNames{ii});
        m(areaNames{ii}) = vol;
        fprintf('%d/%d. %s - %0.3f mm^3\n', ii, length(areaNames), areaNames{ii},vol)
    end

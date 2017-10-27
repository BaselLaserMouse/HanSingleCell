function saveLaminarPlots(lamData,cleanCells,areaFits)
% save all arbor plots


    for ii=1:length(lamData)
        if length(lamData(ii).cellID)<6, continue, end
        arborPlotWithDensity(lamData,cleanCells,[],ii,areaFits{ii});
        drawnow
        fname = [strrep(lamData(ii).areaName, ' ','_'), '_arborPlot'];
        print('-dpng',fullfile('laminarFigs',fname))
    end


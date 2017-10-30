function saveLaminarPlots(lamData,cleanCells,areaFits)
% save all arbor plots


    lamPlotInd = [1,2,3,4,5,9,13];

    for ii=1:length(lamPlotInd)
        %if length(lamData(ii).cellID)<6, continue, end
        ind = lamPlotInd(ii);
        arborPlotWithDensity(lamData,cleanCells,[],ind,areaFits{ind});
        drawnow
        fname = [strrep(lamData(ind).areaName, ' ','_'), '_arborPlot'];
        print('-dpng',fullfile('laminarFigs',fname))
        print('-depsc',fullfile('laminarFigs',fname))

        topViewProjectionMultipleCells(cleanCells, lamData(ind) );
        fname = [strrep(lamData(ind).areaName, ' ','_'), '_PHP'];
        print('-depsc',fullfile('laminarFigs',fname))
    end


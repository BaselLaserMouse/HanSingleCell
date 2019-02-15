function S=overlaySubCortical(H,viewName,subCorticalInd)
    %
    %
    % Overlay sub-cortical areas on to a view:
    % H - handles structure from overlayCellsOnThreeProjections
    % viewName - 'sagittal', 'transverse', 'coronal'
    % subCorticalInd - brain area index values 
    %   e.g. for striatum, amygdala and ACC (multiple layers grouped) : {672,131,[211,935,1015]})

    S=[];
    if isempty(subCorticalInd)
        return
    end

    A=aratools.atlascacher.getCachedAtlas;
    vol = A.atlasVolume;

    switch lower(viewName)
        case 'transverse'
        vol = permute(vol, [3,2,1]);
        case 'sagittal'
            vol = permute(vol, [1,3,2]);
        case 'coronal'
        % nothing
        otherwise
            fprintf('UNKNOWN viewName in buildImages > overlaySubCortical: %s\n', viewName)
            return
    end

    %loop through the sub-cortical area index values and overlay them on the plot
    for ii=1:length(subCorticalInd)
        these_areas = subCorticalInd{ii};
        areaMask = any(vol==these_areas(1),3);

        if length(these_areas)>1
            for kk=2:length(these_areas)
                areaMask = areaMask+any(vol==these_areas(kk),3);
            end
        end
        areaMask = any(areaMask,3);

        viewName = [upper(viewName(1)), lower(viewName(2:end))];

        set(gcf,'CurrentAxes',H.(['axes',viewName])) %Change focus
        hold on
        B=bwboundaries(areaMask);
        for kk=1:length(B)
            S(end+1)=plot(B{kk}(:,2),B{kk}(:,1),'b--');
        end
        hold off
    end



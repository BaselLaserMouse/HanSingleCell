function S=overlaySelectAreas(H,ARAindex)
    % buildImages helper function
    %
    %Only overlay areas that this neuron projects to
    areasToOverlay={};
    if any(ARAindex==211) || any(ARAindex==935) || any(ARAindex==1015) || ...
        any(ARAindex==588) || any(ARAindex==296) || any(ARAindex==772) || ...
        any(ARAindex==810) || any(ARAindex==919)
        %ACC
        areasToOverlay{end+1} = [211,935,1015,588,296,772,810,919];
    end

    if any(ARAindex==131) 
        %lateral amygdalar
        areasToOverlay{end+1} = 131;
    end

    if any(ARAindex==672) 
        %lateral amygdalar
        areasToOverlay{end+1} = 672;
    end

    S=[overlaySubCortical(H,'transverse',areasToOverlay),...
        overlaySubCortical(H,'sagittal',areasToOverlay),...
        overlaySubCortical(H,'coronal',areasToOverlay)];


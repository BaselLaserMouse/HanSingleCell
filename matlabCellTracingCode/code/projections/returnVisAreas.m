function [visAreas,moreAreas] = returnVisAreas(outlinesDir)
    % Return a cell array of paths to (mainly) visual area outlines. You can feed these to 
    % overlayManyCellsOnThreeProjections as the second argument.
    %
    % function [visAreas,moreAreas] = returnVisAreas(outlinesDir)
    % 
    %
    % Inputs
    % outlinesDir is: ~/tvtoucan/Mrsic-Flogel/ReferenceAtlas/ARA_CCFv3/ARA_25_micron_mhd/outlines_smoothed
    %                   by default.
    %
    % Outputs
    % visAreas - the areas close to V1
    % moreAreas - mode distant areas that also recieve projections from V1
    %
    %
    % Distinction between the above two is somewhat arbirtary. Based on what plots well at the moment. 
    %
    %
    % Also See:
    % colorBrains


if nargin<1
    outlinesDir='~/tvtoucan/Mrsic-Flogel/ReferenceAtlas/ARA_CCFv3/ARA_25_micron_mhd/outlines_smoothed';
end


areaNamesV = {'AA.mat', 'LIA.mat', 'PLVA.mat', ...
             'ALVA.mat', 'LVA.mat', 'PMVA.mat', ...
             'RLV.mat', 'AMVA.mat'};

areaNamesM = {'postrhinal.mat',  'TEA.mat', 'retrosplenial.mat', 'ectorhinal.mat', 'entorhinal.mat'};



for ii=1:length(areaNamesV)
    visAreas{ii} = fullfile(outlinesDir,areaNamesV{ii});
end


for ii=1:length(areaNamesM)
    moreAreas{ii} = fullfile(outlinesDir,areaNamesM{ii});
end



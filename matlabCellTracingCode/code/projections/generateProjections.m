function [out,ATLAS] = generateProjections(ATLAS,surfaceDepth)
% Generates projections for the three axes 
%
% function OUT = generateProjections(ATLAS,surfaceDepth)
%
% Purpose
% The structure produced by this function is saved to disk so it can be
% used by the associated plotting functions in this directory. 
% The structure contains the projected atlas along all three axes.
%
% Inputs
% ATLAS - atlas volume
% surfaceDepth - vector of length three defining the surfaceDepth parameter for 
%                for: [top-down,side,front]. See help aratools.projectAtlas.generate
%                default: [4,12,12]
%
% Example
% ATLAS=mhd_read('~/tvtoucan/Mrsic-Flogel/ReferenceAtlas/ARA_CCFv3/ARA_25_micron_mhd/atlas_smooth1_corrected.mhd');
% out = generateProjections(ATLAS);
%
% 
% Also see the following functions from ARAtools:
% aratools.projectAtlas.generate
% aratools.projectAtlas.createBoundaries
% aratools.projectAtlas.highlightAreaPlot
% 
%
% Rob Campbell - Basel 2017


if nargin<2
    surfaceDepth=[4,3,9];
end

%Merge the bulb into one thing
%ATLAS(ATLAS==1016)=507;

%Merge motor areas and surrounding
%ATLAS(ATLAS==320)=656;
%ATLAS(ATLAS==667)=656;

grouping={[656,320,667,68,943,171,484,36,962,204],... %motor and around
          [1009,477,707,588,776],...
          [442,542,73,434,430],... %retroplenial areas
          [1017,1025,944],... %tidy cerebellum
          [1030,113],...
          [1006,670],...
          [1074,905],...
          [1007,313],...
          [735,251],...
          [558,838],...
          [0,1049],...
          [593,821],... %V1
          [805,41],... %higher visual
          [959,1127,755,696],... %auditory
          [981,201],...
          [750,269],...
          [507,1016],... %bulb
          [836,427],... %merge ectorhinal layers
          [918,926],... %Entorhinal area with Entorhinal area, medial part, dorsal zone
          [996,704,120],... % agranular insular
          };

cerebellumAndBrainStem = {'Cerebellum','Pons','Medulla','arbor vitae','medial longitudinal fascicle',...
'fiber tracts','inferior cerebellar peduncle','ventricular systems'};

cleanCerebellum={'Paraflocculus','arbor vitae'};

groupBulb={'Main olfactory bulb'};
%Top
ii=1;
fprintf('Generating view %d\n',ii)

out(ii) = aratools.projectAtlas.generate(ATLAS,ii,'surfaceDepth',surfaceDepth(ii),'removeAreaWithChildren',...
    cleanCerebellum,'groupIndsInProjection',grouping,'dilateSize',4);

%Side
ii=2;
fprintf('Generating view %d\n',ii)
cerebellumBrainStem2=['Brain stem',cerebellumAndBrainStem];

moreGrouping={...
            [354,978],...
            [235,863],... %brainstem
            [101,403],...
            [771,380,658,841,235,78,114,1009,153,661,1069,190,307,83,1123,794,413,354,429],...
            [4,1007],... %IC
            [961,788],...
            [512,728],...
            [918,566,647,1089]};

groupingSide = [grouping,moreGrouping];

out(ii) = aratools.projectAtlas.generate(ATLAS,ii,'surfaceDepth',surfaceDepth(ii),'groupChildren','Cerebellum',...
    'groupIndsInProjection',groupingSide,'dilateSize',2); 



%Rear view
ii=3;

yetMore = {...
        [354,1123,235,403,1098,136,794,1039,711,903],...
        [957,728],...
        [153,101],....
        [918,961,788,566,647,1089,639,260,754,1009,23,413,], ... %repeat of entorhinal
};
groupingBack = [grouping,yetMore];

fprintf('Generating view %d\n',ii)
out(ii) = aratools.projectAtlas.generate(ATLAS,ii,'surfaceDepth',surfaceDepth(ii),'dilateSize',2,'groupIndsInProjection',groupingBack)%,'removeAreaWithChildren',cerebellumBrainStem2);


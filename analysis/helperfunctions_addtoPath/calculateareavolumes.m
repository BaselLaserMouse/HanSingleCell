%% function to measure the volume of brain areas using the ARA.
function [V]=calculateareavolumes(areas)
%load atlas
load ARA.mat

for i=1:length(areas)
    S=getAllenStructureList('childrenOf',areas(i));
    HIT=ismember(ARA,table2array(S(:,1)));
    V(i)=sum(sum(sum(HIT)));
end
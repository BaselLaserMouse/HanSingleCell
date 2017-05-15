%% function to simulate double-retrograde tracing based on a certain projection strength cutoff, that is required for "labeling" with the tracer.
function [pairs,counts,overlap]=simulatedretrograde(matrix,cutoff,labels)
matrix_binary=matrix>=cutoff;
x=1:size(matrix_binary,2);
motifs=nchoosek(1:size(matrix_binary,2),2);

for i=1:size(motifs,1)
    hitone=sum(matrix_binary(:,motifs(i,1))~=0);
    hittwo=sum(matrix_binary(:,motifs(i,2))~=0);
	hitboth=sum(matrix_binary(:,motifs(i,1))==1&matrix_binary(:,motifs(i,2))==1);
    
    counts(i,:)=[hitone,hittwo,hitboth];
    overlap(i,:)=[hitboth/(hitone+hittwo-hitboth),hitboth/hitone,hitboth/hittwo]*100;
    pairs(i,:)=labels(motifs(i,:));
end

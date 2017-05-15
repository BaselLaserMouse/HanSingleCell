%% function to calculate and plot a set of first order stats
function [ded]=firstorderstats_ratios(matrix,ratio,cutoff,labels,area_volumes)

matrix_tmp=matrix(sum(matrix>=cutoff,2)~=0,:);
matrix_f=matrix_tmp(:,sum(matrix_tmp>=cutoff,1)~=0);
labels=labels(sum(matrix_tmp>=cutoff,1)~=0);

matrix_f(matrix_f<cutoff)=0;
% matrix_f=matrix_f./repmat(area_volumes,size(matrix_f,1),1);
matrix_fn=matrix_f./repmat(max(matrix_f,[],2),1,size(matrix_f,2));


targets=matrix_fn>ratio;
targets_num=sum(targets,2);
input_num=sum(targets,1);

%overall distribution of target numbers.
h_overalldistribution=hist(targets_num',1:size(matrix_f,2));

%fraction of dedicated input to each area
ded=h_overalldistribution(1)/sum(h_overalldistribution);

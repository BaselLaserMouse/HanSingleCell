%% function to compare the distances between a reference set of neurons and a testcase
function [out,out_ref,out_test,out_random]=sampleandcomparepoints(reference,test,repeats,distance)
out=zeros(repeats,3);
out_test=zeros(repeats,size(test,1));
out_ref=zeros(repeats,size(test,1));
out_random=zeros(repeats,size(test,1));


for R=1:repeats
%sample reference
r=randsample(size(reference,1),size(test,1));
sample=reference(r,:);
x=ones(size(reference,1),1);
x(r)=0;
remains=reference(logical(x),:);

%generate random neurons by sampling from a uniform distribution
random_tmp=rand(size(remains));
randomneurons=random_tmp./repmat(sum(random_tmp,2),1,size(random_tmp,2));
% alternatively generate random neurons by reshuffling area labels for
% every neuron
% r2=randsample(size(reference,1),size(test,1));
% sample2=reference(r2,:);
% 
% randomneurons=sample2(randi(size(sample2,2),size(sample2)));


% calculate distance of referencesample to closest remains datapoit
[D_s_tmp]=pdist2(sample,remains,distance);
[D_s,I_s]=min(D_s_tmp,[],2);


% calculate distance of testsample to closest remains datapoint
[D_t_tmp]=pdist2(test,remains,distance);
[D_t,I_t]=min(D_t_tmp,[],2);


% calculate distance of test to closest random neuron
[D_r_tmp]=pdist2(test,randomneurons,distance);
[D_r,I_r]=min(D_r_tmp,[],2);



%prepare output
out(R,:)=[mean(D_s),mean(D_t),mean(D_r)];
out_test(R,:)=D_t;
out_ref(R,:)=D_s;
out_random(R,:)=D_r;
end



%% function to collapse projections of individual neurons to adjacent areas into one. very conservative, as chaining of areas is counted as one
% e.g. projections to PM,AM,A and RL will be counted as 1.
function []=removeneighbours(matrix,cutoff,area_id)

matrix_tmp=matrix(sum(matrix>=cutoff,2)~=0,:);
matrix_f=matrix_tmp(:,sum(matrix_tmp>=cutoff,1)~=0);
area_id=area_id(sum(matrix_tmp>=cutoff,1)~=0);

targets=matrix_f>=cutoff;
targets_num=sum(targets,2);

%find all the neighbouring areas.
load leavesbyArea.mat
neighbours=zeros(size(matrix_f,2));
for i=1:size(matrix_f,2);
     [~,U]=X.getBorderPixelsForArea(area_id(i),true);
     [d,loc]=ismember(U,area_id);
     areas(i).neighbours=loc;
     neighbours(loc(loc~=0),i)=1;
end
figure;


% find all connected components and collapse into one.
for i=1:size(matrix_f,1)
    neuron_neighbours=neighbours;
        neuron_neighbours(targets(i,:)==0,:)=0;
        neuron_neighbours(:,targets(i,:)==0)=0;
        neuron_neighbours=neuron_neighbours(sum(neuron_neighbours,2)~=0,sum(neuron_neighbours,1)~=0);
        if size(neuron_neighbours)~=[0 0]
        [S(i),G]=graphconncomp(sparse(neuron_neighbours),'Directed','false');
        G_size(i)=length(G);
        else
            S(i)=0
            G_size(i)=0;
        end
end
targets_adjust=targets_num-G_size'+S';
pie(hist(targets_adjust',1:size(matrix_f,2)));
figure;
pie(hist(targets_adjust(targets_num>1)',1:size(matrix_f,2)));


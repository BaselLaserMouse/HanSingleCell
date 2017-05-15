%% function to calculate and plot a set of first order stats
function []=firstorderstats(matrix,cutoff,labels)

%prepare the projection matrix
matrix_tmp=matrix(sum(matrix>=cutoff,2)~=0,:); %remove cells that do not cross projection cutoff for any target area
matrix_f=matrix_tmp(:,sum(matrix_tmp>=cutoff,1)~=0); %remove target areas that receive no input
labels=labels(sum(matrix_tmp>=cutoff,1)~=0); %trim the area labels accordingly
targets=matrix_f>=cutoff; %binarize the projection matrix

%calculate raw stats
targets_num=sum(targets,2); %calculate the number of targets per neuron
input_num=sum(targets,1); %calculate the total number of neurons projecting to each target area

% plot the overall distribution of target numbers.
h_overalldistribution=hist(targets_num',1:size(matrix_f,2));

figure;% in the shape of a bar graph
bar(1:size(matrix_f,2),h_overalldistribution,0.5);
box off
xlabel('Number of projection targets')
ylabel('Number of cells')

figure; % or in the shape of a pie chart
explode=[1, zeros(1,size(matrix_f,2)-1)];
pie(h_overalldistribution,explode)
title('Overall distribution of target numbers')
box off


%plot the fraction of dedicated input to each area

figure; %overall inneravation pattern
innervation=input_num./size(matrix_f,1);
[innervation_s,ix]=sort(innervation,'descend');
bar(1:size(matrix_f,2),innervation_s,0.5);
ax = gca;
ax.XTick=[1:size(matrix_f,2)];
ax.XTickLabel=labels(ix);
rotateXLabels( gca(), 90 )
xlabel('Target area')
ylabel('Fraction of input cells')
title('Fraction of input to each area')
box off


frac_dedicated=[];
for areas=1:size(matrix_f,2)
    total=input_num(areas);
    dedicated=sum(targets(:,areas)==1 & targets_num==1);
    frac_dedicated(areas)=dedicated./total;
end

figure; %fraction of dedicated input
bar(1:size(matrix_f,2),frac_dedicated(ix),0.5)
ax = gca;
ax.XTick=[1:size(matrix_f,2)];
ax.XTickLabel=labels(ix);
rotateXLabels( gca(), 90 )
xlabel('Target area')
ylabel('Fraction of dedicated input cells')
title('Fraction of dedicated input to each area')
box off

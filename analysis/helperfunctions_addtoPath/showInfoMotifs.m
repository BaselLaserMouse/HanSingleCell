%% function to show all the interesting numbers for identified statistically significant motifs.
function []=showInfoMotifs(matrix_raw,origin,cutoff,significantmotifs,directionofregulation,observed,expected,pvalue,labels,plotall)
if plotall~=1
figure;
end
load brainAreaProjections.mat %load dataobject that allows to draw ARA outlines

cells=matrix_raw>=cutoff; %binarize projections
cells_norm=matrix_raw./repmat(max(matrix_raw,[],2),1,size(matrix_raw,2)); %normalize to the max

for i=1:size(significantmotifs,1);
    if plotall==1
        figure;
    else
    clf
    end
    
    d=ismember(cells,significantmotifs(i,:),'rows');
    h=hist(origin(d),4:7);
    freq=h./hist(origin,4:7);
    
    %plot
    subplot(2,2,1)
    make_nice_MAPseqplots3(projections,significantmotifs(i,:),8,2);
    subplot(2,2,2);
    ratio=observed(i)/expected(i);
    text(0.1,0.9,['regulation ',num2str(directionofregulation(i))]); axis off
    text(0.1,0.7,['log2 ratio ',num2str(log2(ratio))]); 
    text(0.1,0.5,['pvalue ',num2str(pvalue(i))]); 
    text(0.1,0.3,['observed ',num2str(observed(i)),'; expected ',num2str(expected(i))]); 

    %plot the projection strengths of all neurons projecting like this,
    %color coded by mouse of origin
    subplot(2,2,3);
    cells_norm_tmp=cells_norm(d,:);
    origin_tmp=origin(d);
    mousecolor=[       0    0.4470    0.7410;0.8500    0.3250    0.0980;...
    0.9290    0.6940    0.1250;...
    0.4940    0.1840    0.5560];
    for m=unique(origin)
        plot(cells_norm_tmp(origin_tmp==m,:)','o-','Color',mousecolor(m-3,:));
        hold on
    end
    ylabel('Projection strength')
    ax = gca;
    ax.XTick=[1:length(labels)];
    ax.XTickLabel=labels;
    box off
    
    subplot(2,2,4);
    cells_tmp=matrix_raw(d,:);
    origin_tmp=origin(d);
    mousecolor=[       0    0.4470    0.7410;0.8500    0.3250    0.0980;...
    0.9290    0.6940    0.1250;...
    0.4940    0.1840    0.5560];
    for m=unique(origin)
        plot(cells_tmp(origin_tmp==m,:)','o-','Color',mousecolor(m-3,:));
        hold on
    end
    ylabel('Barcode count')
    ax = gca;
    ax.XTick=[1:length(labels)];
    ax.XTickLabel=labels;
    box off
    
    
    if plotall~=1
    pause();
    end
end
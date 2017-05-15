function dostatsforMAPseqneurons
%% 
% Justus Kebschull

%% load V1 projection matrix
load matrix_norm.mat
matrix_raw=matrix_norm(:,2:end);

cutoff=5; %binarization threshold; at least 5 barcode counts per area to be considered a projection
labels={'Li','LM','AL','PM','AM','RL'};


% load area volumes
load V6areavolumes.mat

%% make a clustergram showing the projection matrix
%load the white-to-orange colormap
load heatmap_colormap.mat
% normalize projection matrix by each neuron's maximum projection
matrix_clust=matrix_norm./repmat(max(matrix_norm,[],2),1,size(matrix_norm,2));
labels_clust={'OB',labels{:}};

clustergram(matrix_clust,...
    'Symmetric',0,'Standardize',3,'Cluster',1,'Colormap',cmap2,'ColumnLabels',labels_clust,...
    'RowPDist','cosine','ColumnPDist','cosine');

%% calculate first order statistics
firstorderstats(matrix_raw,cutoff,labels);

%% calculate the fraction of broadcasting cells enforcing a minimum ratio between maxium and other projections
counter=1;
for ratio=0:0.1:1;
dedicated_tmp=firstorderstats_ratios(matrix_raw,ratio,cutoff,labels,V);
dedicated(counter)=dedicated_tmp(1);
counter=counter+1;
close all;
end

figure;
plot(0:0.1:1,1-dedicated);
box off
xlabel('Minimum ratio of secondary projection to maximum')
ylabel('Fraction of broadcasting cells');


%% plot conditional probabilities
[conditionalp,rowlabels,columnlables,matrix_binary]=conditionalprobabilities(matrix_raw,5,cutoff,labels,0);
sum(matrix_binary,1)


conditionalp2=conditionalp-2*eye(length(labels));
figure;imagesc(conditionalp2');
    ax = gca;
    ax.YTick=[1:length(columnlables)];
    ax.YTickLabel=columnlables;
    ax.XTick=[1:length(rowlabels)];
    ax.XTickLabel=rowlabels;
    rotateXLabels(ax,45)
    ylabel('P(B|A)')
    xlabel('Area A')
    caxis([0 0.8]);




%% look at polyfrucations
[motif,pvalue_all,significance,overorunder,observed_all,expected_all,p_expected,neuronnum]=higherorder(matrix_raw,cutoff,0.05)
observed=observed_all(significance==1);
expected=expected_all(significance==1);
directionofregulation=overorunder(significance==1);
pvalue=pvalue_all(significance==1);
significantmotifs=motif(significance==1,:);

figure;
load heatmap_colormap.mat
colormap(cmap2)
imagesc(motif')

%% check the mouse of origin for each of the overrepresented motifs
showInfoMotifs(matrix_raw,origin,cutoff,significantmotifs,directionofregulation,...
    observed,expected,pvalue,labels,1);




%% plot scatter plot of LM-AL bifurcation
matrix_binary=matrix_raw>=cutoff;
matrix_norm=matrix_raw./repmat(max(matrix_raw,[],2),1,6);
LMAL=ismember(matrix_binary,[0 1 1 0 0 0],'rows');
LMAL=matrix_binary(:,2)==1&matrix_binary(:,3)==1;
figure;
c=parula(4);

for o=4:7;
    subplot(2,2,o-3);
    scatter(matrix_raw(LMAL&origin'==o,2),matrix_raw(LMAL&origin'==o,3),25);
    title(['mouse ',int2str(o)]);
    xlabel('LM counts')
    ylabel('AL counts')
end

%% how much overlap would we see if we had done double retrograde tracing?
cutoff=5;
[pair,counts,overlap]=simulatedretrograde(matrix_raw,cutoff,labels);    

O=table(pair,counts,overlap)
    



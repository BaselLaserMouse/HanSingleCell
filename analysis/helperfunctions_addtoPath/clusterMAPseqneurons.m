function clusterMAPseqneurons

%% this script will cluster MAPseq neurons 
% we will attempt to cluster the MAPseq neurons using kmeans clustering and
% will give examples of each cluster of neurons using the raster view of
% the visual cortex space.
%
% Justus Kebschull



load matrix_norm.mat

%% normalize data
matrix_tmp=matrix_norm(:,[2:7]);
matrix=matrix_tmp./repmat(max(matrix_tmp,[],2),1,size(matrix_tmp,2));


%% check cluster stability
clusterstability_kmeans=evalclusters(matrix,'kmeans','gap','Klist',[1:20],'Distance','cosine');
figure;plot(clusterstability_kmeans.CriterionValues);
xlabel('number of clusters')
ylabel('GAP criterion')
box off
xlim([1 20])

clusterstability_kmeans2=evalclusters(matrix,'kmeans','silhouette','Klist',[1:20],'Distance','cosine');
figure;plot(clusterstability_kmeans2.CriterionValues);
xlabel('number of clusters')
ylabel('silhouettee criterion')
box off


%% cluster MAPseq data
k=8; %chosen cluster number

[T,C,sumd,D]=kmeans(matrix,k,'Distance','cosine');

%% make a nice plot

maketiledexampleplot(k,T,C,7,matrix,1);

%% save plots to individual files
E=8;
load brainAreaProjections.mat

for cluster=1:k;
    neurons=matrix(T==cluster,:);
    figure;
    for example=1:E/2
        subplot(2,2,example)
        make_nice_MAPseqplots3(projections,neurons(example,:),7,3);
    end
        print(['examples\',int2str(cluster),'_',int2str(1)],'-dpdf')
    figure;
    for example=E/2+1:E
        subplot(2,2,example-E/2)
        make_nice_MAPseqplots3(projections,neurons(example,:),7,3);
    end
        print(['examples\',int2str(cluster),'_',int2str(2)],'-dpdf')
    close all
end

figure;

for cluster=1:8;
    neurons=matrix(T==cluster,:);
    subplot(4,2,cluster)
        make_nice_MAPseqplots3(projections,C(cluster,:),7,3);
        %close all
    
end
print(['examples\cluster_cluster2'],'-dpdf')


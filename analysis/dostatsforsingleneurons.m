function dostatsforsingleneurons
%% analyse single neuron dataset
% Justus Kebschull


%% load data

load OUT_2updatedlength.mat %load dataset with an exclusion zone of 50um on each side of the border
% load OUT_2updatedlength_all.mat %load dataset with an exclusion zone of % 50um on each side of the border - if looking at non-v1 cells


load goodcells.mat% load the names of all the high confidence fills that will remain in analysis
load nonV1.mat %load the names of all the good fills of cells that do not reside in V1

cutoff=200; %set minimum length cutoff for an area to be considered a target area. unit here is 5um, so 200*5=1000um
labels=OUT.areaNamesInSamples; %extract area labels

%remove V1 from the matrix
v1=strmatch('Primary visual area',labels);
x=ones(length(labels),1);
x(v1)=0;
labels=labels(logical(x))






%keep only layer 2/3 cells
matrix_tmp=OUT.dataMat(logical(x),OUT.somaLocation==821 | OUT.somaLocation==593 | OUT.somaLocation==721)';
area_ids=OUT.allAreas(logical(x),:);
OUT.cellIDs=OUT.cellIDs(OUT.somaLocation==821 | OUT.somaLocation==593 | OUT.somaLocation==721);

% matrix_tmp=OUT.dataMat(logical(x),:)';
% area_ids=OUT.allAreas(logical(x),:);

%remove striatum retrograde cells
retrograde=goodbackcells;
% retrograde=nonv1cells; %alternatively use these variables to make plots
% of the non-V1 cells
counter=1;
m=[];
for i=1:length(retrograde)
    try
    m(counter)=strmatch(retrograde(i),OUT.cellIDs);
    counter=counter+1;
    catch
    end
end
x=ones(1,size(matrix_tmp,1));
x(m)=0;
matrix_tmp2=matrix_tmp(logical(x),:);
matrix_back=matrix_tmp(~logical(x),:); % projection matrix for all backlabeled cells

cellIDs=OUT.cellIDs(logical(x));

%remove abrupt terminations, i.e. everything that is not a "good" fill.

counter=1;
m=[];
for i=1:length(goodcells)
    try
    m(counter)=strmatch(goodcells(i),cellIDs);
    counter=counter+1;
    catch
    end
end
x=ones(1,size(matrix_tmp2,1));
x(m)=0;
matrix_tmp3=matrix_tmp2(~logical(x),:);


matrix_raw=matrix_tmp3(sum(matrix_tmp3,2)~=0,:); % projection matrix for all randomly labeled V1 cells

%% restrict to six hihger visual areas (LM, AL, AM, PM, RL, LI) -- chose either this block OR the next block, not both to run this script.
%reduce projection matrix to only contain "visual areas"
a={brainAreaNames.visualAreas{:}};
m=[];
counter=1
for i=1:length(brainAreaNames.visualAreas)
    try
    m(counter)=strmatch(a(i),labels);
    counter=counter+1;
    catch
    end
end

matrix_raw=matrix_raw(:,m);
labels=labels(m);
area_ids=area_ids(m,:);

%further reduce the projection matrix to contain only the six relevant
%higher visual areas.
match=[12 4 2 3 1 10];
matrix_raw=matrix_raw(:,match);
labels=labels(match);

matrix_raw=matrix_raw(sum(matrix_raw,2)~=0,:); %final projection matrix


%% restrict analysis to potential target areas - no outofBrain, no whitematter. -- chose either this block OR the preceeding block, not both to run this script.

%find all the non-grey matter areas
[ind,whitematter_names]=aratools.utils.whiteMatterInds;
a={whitematter_names{:},'Out of brain','lateral ventricle','corpus callosum','fiber tracts','dorsal hippocampal commissure','alveus'};
counter=1
m=[];
for i=1:length(a)
    try
    m(counter)=strmatch(a(i),labels);
    counter=counter+1;
    catch
    end
end
x=ones(length(labels),1);
x(m)=0;

%remove these areas from all projection matrices and from the area labels.
matrix_raw=matrix_raw(:,logical(x));
matrix_back=matrix_back(:,logical(x));
labels=labels(logical(x));
area_ids=area_ids(logical(x),:);
matrix_raw=matrix_raw(sum(matrix_raw,2)~=0,:);

%% calculate area volumes by counting up the voxels of the atlas
[V]=calculateareavolumes(labels);


%% calculate first order statistics

firstorderstats(matrix_raw,cutoff,labels); %for the randomly selected V1 cells
% firstorderstats(matrix_back,cutoff,labels); % for the backlabeled cells

%% calculate the fraction of broadcasting cells enforcing a minimum ratio between maxium and other projections
counter=1;
for ratio=0:0.05:1;
dedicated_tmp=firstorderstats_ratios(matrix_raw,ratio,cutoff,labels,V);
dedicated(counter)=dedicated_tmp(1);
counter=counter+1;
close all;
end

plot(0:0.05:1,1-dedicated);
box off
xlabel('Minimum ratio of secondary projection to maximum')
ylabel('Fraction of broadcasting cells');

%% make a nice clustergram
% enforce order:
a={brainAreaNames.visualAreas{:}};

[m,loc]=ismember(labels,a);
[s,ix]=sort(loc);
relabel=[ix(sum(s==0)+1:end);ix(1:sum(s==0))];
labels=labels(relabel);
matrix_raw=matrix_raw(:,relabel);
matrix_back=matrix_back(:,relabel);

%normalize projection matices
matrix_clust_tmp=matrix_raw./repmat(max(matrix_raw,[],2),1,size(matrix_raw,2));%normalize by max
matrix_back_clust_tmp=matrix_back./repmat(max(matrix_back,[],2),1,size(matrix_back,2));%normalize by max

%remove areas that do not receive cross-threshold input from at least one cell 
areastokeep=sum([matrix_raw;matrix_back]>=cutoff,1)~=0;
matrix_clust=matrix_clust_tmp(:,areastokeep);

labels_clust=labels(areastokeep);
matrix_back_clust=matrix_back_clust_tmp(:,areastokeep);



load heatmap_colormap.mat
%plot a clustered projection heatmap for the randomly selected V1 cells
cg=clustergram(matrix_clust(sum(matrix_raw>=cutoff,2)~=0,:),...
    'Symmetric',0,'Standardize',3,'Cluster',1,'Colormap',cmap2,'ColumnLabels',labels_clust,...
    'RowPDist','cosine','ColumnPDist','cosine');
%plot a clustered projection heatmap for the backlabeled cells
cg=clustergram(matrix_back_clust(sum(matrix_back>=cutoff,2)~=0,:),...
    'Symmetric',0,'Standardize',3,'Cluster',1,'Colormap',cmap2,'ColumnLabels',labels_clust,...
    'RowPDist','cosine','ColumnPDist','cosine');

%% plot figure showing resilence to parameter changes

% buffer zone size; fixed 200 length cutoff
dedicated_buffer=[20 20 23 25 29 31 45 55 76];
bufferzone=[0 25 50 75 100 125 150 200 250];

figure;
plot(bufferzone,100-dedicated_buffer,'o-')
xlabel('size of bufferzone (um)')
ylabel('% broadcasting cells')
ylim([0 100])
box off

% cutoff; buffer fixed at 50um
cut=[1 50:50:500];
dedicated_cutoff=[10 10 19 20 23 23 21 25 25 39 43];

figure;
plot(cut*5/1000,100-dedicated_cutoff,'o-');
xlabel('minimum legth per traget (mm)')
ylabel('% broadcasting cells')
ylim([0 100])

box off


%% look at fraction of broadcasting cells that project to non-neighbouring areas.
removeneighbours(matrix_raw,cutoff,area_ids(:,1));
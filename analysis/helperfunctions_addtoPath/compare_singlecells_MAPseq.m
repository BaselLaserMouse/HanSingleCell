function compare_singlecells_MAPseq
% Justus Kebschull
%% compare single celll dataset to MAPseq data



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

%% restrict to six higher visual areas (LM, AL, AM, PM, RL, LI) -- chose either this block OR the next block, not both to run this script.
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

singleneurons_tmp2=matrix_raw(sum(matrix_raw,2)~=0,:); %final projection matrix

singleneurons=singleneurons_tmp2./repmat(max(singleneurons_tmp2,[],2),1,size(singleneurons_tmp2,2));


%% load MAPseq data 
load matrix_norm.mat
cutoff=5;
matrix_tmp=matrix_norm(sum(matrix_norm(:,[2:end]),2)~=0,[2:end]);
matrix=matrix_tmp./repmat(max(matrix_tmp,[],2),1,size(matrix_tmp,2));


%% do analysis independent of centroids


[out,out_ref,out_test,out_random]=sampleandcomparepoints(matrix,singleneurons,1000,'cosine');

h1=hist(out_ref(:),0:0.01:0.3);
h2=hist(out_test(:),0:0.01:0.3);
h3=hist(out_random(:),0:0.01:0.3);



figure;plot(0:0.01:0.3,h1/sum(h1));hold on;plot(0:0.01:0.3,h2/sum(h2));plot(0:0.01:0.3,h3/sum(h3));
xlabel('Distance')
ylabel('Probability')
box off
xlim([0 0.3])
% 
[h,p]=kstest2(h1,h2)
[h,p]=kstest2(h1,h3)
[h,p]=kstest2(h2,h3)

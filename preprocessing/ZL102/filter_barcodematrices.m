function filter_barcodematrices
% Justus Kebschull
%%



%% load data

load brain6.mat
load brain7.mat
load spikes.mat

%% filter data;
lowercutoff=10;
uppercutoff=1000;

matrix6f=matrix6(max(matrix6(:,[1:5,8:size(matrix6,2)]),[],2)>lowercutoff,:);
barcodes6f=barcodes6(max(matrix6(:,[1:5,8:size(matrix6,2)]),[],2)>lowercutoff,:);
matrix7f=matrix7(max(matrix7(:,[1:5,8:size(matrix7,2)]),[],2)>lowercutoff,:);
barcodes7f=barcodes7(max(matrix7(:,[1:5,8:size(matrix7,2)]),[],2)>lowercutoff,:);

%normalize the matrices to spike ins and to sum to one
x=[];for i=1:36;x(i)=length(spikes(i).counts2u);end
matrix6n_tmp=matrix6f./repmat(x(1:17),size(matrix6f,1),1);
   
matrix7n_tmp=matrix7f./repmat(x(18:34),size(matrix7f,1),1);
    
%% check for things with cell body in V1 for mouse6
% require for good cells to have much more projection to V1 than elsewhere.
minratio=0.1; %minum ratio of second highest, and highest expression level
minpeak=300; %minimum molecule count for a cell body
[cellbodies6,matrix6ff,dump,keep6]=findcellbodies(matrix6n_tmp,minratio,barcodes6f);
rawcounts6_f=matrix6f(keep6,:);

[m,loc]=max(matrix6ff,[],2);
matrix6_v1upper_tmp=matrix6ff((loc==6 & rawcounts6_f(:,6)>minpeak),[1:5,8:17]);
matrix6_v1upper=matrix6_v1upper_tmp./repmat(sum(matrix6_v1upper_tmp,2),1,15);
matrix6_v1lower_tmp=matrix6ff((loc==7 & rawcounts6_f(:,7)>minpeak),[1:5,8:17]);
matrix6_v1lower=matrix6_v1lower_tmp./repmat(sum(matrix6_v1lower_tmp,2),1,15);

figure;hist(loc,0:18) %plot the position of detected cell bodies

matrix6_v1upper_raw=rawcounts6_f((loc==6 & rawcounts6_f(:,6)>minpeak),[1:5,8:17]);
matrix6_v1lower_raw=rawcounts6_f((loc==7 & rawcounts6_f(:,7)>minpeak),[1:5,8:17]);


labels6={'TeA/Ent post','TeA/Ent ant','OB','SC','Li','LM','AL','PM','AM','RL','LP','LGN',...
    'Str','LeftV','MM'}
clustergram(matrix6_v1upper,...
    'Symmetric',0,'Standardize',3,'Cluster',1,'Colormap','parula','ColumnLabels',labels6);
clustergram(matrix6_v1lower,...
    'Symmetric',0,'Standardize',3,'Cluster',1,'Colormap','parula','ColumnLabels',labels6);



%% check for things with cell body in V1 for mouse7
% require for good cells to have much more projection to V1 than elsewhere.
minratio=0.1;
minpeak=300;

[cellbodies7,matrix7ff,dump,keep7]=findcellbodies(matrix7n_tmp,minratio,barcodes7f);
rawcounts7_f=matrix7f(keep7,:);

[m,loc]=max(matrix7ff,[],2);
matrix7_v1upper_tmp=matrix7ff((loc==6 & rawcounts7_f(:,6)>minpeak),[1:5,8:17]);
matrix7_v1upper=matrix7_v1upper_tmp./repmat(sum(matrix7_v1upper_tmp,2),1,15);
matrix7_v1lower_tmp=matrix7ff((loc==7 & rawcounts7_f(:,7)>minpeak),[1:5,8:17]);
matrix7_v1lower=matrix7_v1lower_tmp./repmat(sum(matrix7_v1lower_tmp,2),1,15);
figure;hist(loc,0:17) %plot the position of detected cell bodies


matrix7_v1upper_raw=rawcounts7_f((loc==6 & rawcounts7_f(:,6)>minpeak),[1:5,8:17]);
matrix7_v1lower_raw=rawcounts7_f((loc==7 & rawcounts7_f(:,7)>minpeak),[1:5,8:17]);


labels7={'TeA/Ect_post','TeA/Ect ant','OB','Sc','Li','LM','RL','PM','AM','AL',...
    'LP','LGN','Str','LeftV','MM'}
    
clustergram(matrix7_v1upper,...
    'Symmetric',0,'Standardize',3,'Cluster',1,'Colormap','parula','ColumnLabels',labels7);
clustergram(matrix7_v1lower,...
    'Symmetric',0,'Standardize',3,'Cluster',1,'Colormap','parula','ColumnLabels',labels7);


%% save data
save('processedmatrix2.mat','matrix6_v1lower','matrix6_v1upper','matrix7_v1lower','matrix7_v1upper');
save('processedmatrix2_raw.mat','matrix6_v1lower_raw','matrix6_v1upper_raw','matrix7_v1lower_raw','matrix7_v1upper_raw');
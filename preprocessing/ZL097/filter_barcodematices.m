function filter_barcodematrices
% Justus Kebschull

%% load data

load brain4.mat
load brain5.mat
load spikes.mat

%% filter data;
lowercutoff=10;
uppercutoff=1000;

matrix4f=matrix4(max(matrix4(:,[1:6,9:size(matrix4,2)]),[],2)>lowercutoff,:);
barcodes4f=barcodes4(max(matrix4(:,[1:6,9:size(matrix4,2)]),[],2)>lowercutoff,:);
matrix5f=matrix5(max(matrix5(:,[1:4,7:size(matrix5,2)]),[],2)>lowercutoff,:);
barcodes5f=barcodes5(max(matrix5(:,[1:4,7:size(matrix5,2)]),[],2)>lowercutoff,:);

%normalize the matrices to spike ins and to sum to one
x=[];for i=1:36;x(i)=length(spikes(i).counts2u);end
matrix4n_tmp=matrix4f./repmat(x(1:17),size(matrix4f,1),1);
   
matrix5n_tmp=matrix5f./repmat(x(18:33),size(matrix5f,1),1);
    
%% check for things with cell body in V1 for mouse4
% require for good cells to have much more projection to V1 than elsewhere.
minratio=0.1; %minum ratio of second highest, and highest expression level
minpeak=300; %minimum molecule count for a cell body

[cellbodies4,matrix4ff,dump,keep4]=findcellbodies(matrix4n_tmp,minratio,barcodes4f);
rawcounts4_f=matrix4f(keep4,:);



[m,loc]=max(matrix4ff,[],2);
matrix4_v1upper_tmp=matrix4ff(loc==7 & rawcounts4_f(:,7)>minpeak,[1:6,9:17]);
matrix4_v1upper=matrix4_v1upper_tmp./repmat(sum(matrix4_v1upper_tmp,2),1,15);
matrix4_v1lower_tmp=matrix4ff(loc==8 & rawcounts4_f(:,8)>minpeak,[1:6,9:17]);
matrix4_v1lower=matrix4_v1lower_tmp./repmat(sum(matrix4_v1lower_tmp,2),1,15);

figure;hist(loc,0:18) %plot the position of detected cell bodies


matrix4_v1upper_raw=rawcounts4_f((loc==7 & rawcounts4_f(:,7)>minpeak),[1:6,9:17]);
matrix4_v1lower_raw=rawcounts4_f((loc==8 & rawcounts4_f(:,8)>minpeak),[1:6,9:17]);




labels4={'TeA_r','TeA_m','TeA/Ect_c','OB','SC','Li','LM','AL','PM','AM','RL','LP','LGN',...
    'Str','LeftV'}
clustergram(matrix4_v1upper,...
    'Symmetric',0,'Standardize',3,'Cluster',1,'Colormap','parula','ColumnLabels',labels4);
clustergram(matrix4_v1lower,...
    'Symmetric',0,'Standardize',3,'Cluster',1,'Colormap','parula','ColumnLabels',labels4);


%% check for things with cell body in V1 for mouse5
% require for good cells to have much more projection to V1 than elsewhere.
minratio=0.1; %minum ratio of second highest, and highest expression level
minpeak=300; %minimum molecule count for a cell body

[cellbodies5,matrix5ff,dump,keep5]=findcellbodies(matrix5n_tmp,minratio,barcodes5f);
rawcounts5_f=matrix5f(keep5,:);



[m,loc]=max(matrix5ff,[],2);
matrix5_v1upper_tmp=matrix5ff(loc==5 & rawcounts5_f(:,5)>minpeak,[1:4,7:16]);
matrix5_v1upper=matrix5_v1upper_tmp./repmat(sum(matrix5_v1upper_tmp,2),1,14);
matrix5_v1lower_tmp=matrix5ff(loc==6 & rawcounts5_f(:,6)>minpeak,[1:4,7:16]);
matrix5_v1lower=matrix5_v1lower_tmp./repmat(sum(matrix5_v1lower_tmp,2),1,14);
figure;hist(loc,0:17) %plot the position of detected cell bodies


matrix5_v1upper_raw=rawcounts5_f((loc==5 & rawcounts5_f(:,5)>minpeak),[1:4,7:16]);
matrix5_v1lower_raw=rawcounts5_f((loc==6 & rawcounts5_f(:,6)>minpeak),[1:4,7:16]);



labels5={'TeA/Ect_post','OB','Sc','Li','LM','AL','PM','AM','RL',...
    'LP','LGN','Str','LeftV','TeA/Ect_ant'}
    
clustergram(matrix5_v1upper,...
    'Symmetric',0,'Standardize',3,'Cluster',1,'Colormap','parula','ColumnLabels',labels5);
clustergram(matrix5_v1lower,...
    'Symmetric',0,'Standardize',3,'Cluster',1,'Colormap','parula','ColumnLabels',labels5);

%% save data
save('processedmatrix1.mat','matrix4_v1lower','matrix4_v1upper','matrix5_v1lower','matrix5_v1upper');
save('processedmatrix1_raw.mat','matrix4_v1lower_raw','matrix4_v1upper_raw','matrix5_v1lower_raw','matrix5_v1upper_raw');

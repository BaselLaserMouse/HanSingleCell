function makeconsensusmatrix
% Justus Kebschull

%% load all processed matrixes
load ZL097/processedmatrix1.mat
load ZL102/processedmatrix2.mat

load ZL097/processedmatrix1_raw.mat
load ZL102/processedmatrix2_raw.mat

%% load spike ins
load ZL102/spikes.mat
spikes102=spikes;
x102=[];for i=1:length(spikes102);x102(i)=length(spikes102(i).counts2u);end
load ZL097/spikes.mat
spikes097=spikes;
x097=[];for i=1:length(spikes097);x097(i)=length(spikes097(i).counts2u);end
clear spikes;
%% sort columns for consistent target areas and merge matrix
m4=[4 5 6 7 8 9 10 11 12 13 14 15]; m4ect=[1 2 3];
m5=[2 3 4 5 6 7 8 9 10 11 12 13]; m5ect=[1 14];
m6=[3 4 5 6 7 8 9 10 11 12 13 14]; m6ect=[1 2];
m7=[3 4 5 6 10 8 9 7 11 12 13 14]; m7ect=[1 2];

matrix4u=[matrix4_v1upper(:,m4) sum(matrix4_v1upper(:,m4ect),2)];
matrix5u=[matrix5_v1upper(:,m5) sum(matrix5_v1upper(:,m5ect),2)];
matrix6u=[matrix6_v1upper(:,m6) sum(matrix6_v1upper(:,m6ect),2)];
matrix7u=[matrix7_v1upper(:,m7) sum(matrix7_v1upper(:,m7ect),2)];


matrix4u_raw=[matrix4_v1upper_raw(:,m4) sum(matrix4_v1upper_raw(:,m4ect),2)];
matrix5u_raw=[matrix5_v1upper_raw(:,m5) sum(matrix5_v1upper_raw(:,m5ect),2)];
matrix6u_raw=[matrix6_v1upper_raw(:,m6) sum(matrix6_v1upper_raw(:,m6ect),2)];
matrix7u_raw=[matrix7_v1upper_raw(:,m7) sum(matrix7_v1upper_raw(:,m7ect),2)];


matrix4l=[matrix4_v1lower(:,m4) sum(matrix4_v1lower(:,m4ect),2)];
matrix5l=[matrix5_v1lower(:,m5) sum(matrix5_v1lower(:,m5ect),2)];
matrix6l=[matrix6_v1lower(:,m6) sum(matrix6_v1lower(:,m6ect),2)];
matrix7l=[matrix7_v1lower(:,m7) sum(matrix7_v1lower(:,m7ect),2)];

matrix4l_raw=[matrix4_v1lower_raw(:,m4) sum(matrix4_v1lower_raw(:,m4ect),2)];
matrix5l_raw=[matrix5_v1lower_raw(:,m5) sum(matrix5_v1lower_raw(:,m5ect),2)];
matrix6l_raw=[matrix6_v1lower_raw(:,m6) sum(matrix6_v1lower_raw(:,m6ect),2)];
matrix7l_raw=[matrix7_v1lower_raw(:,m7) sum(matrix7_v1lower_raw(:,m7ect),2)];


matrix=[matrix4u;matrix5u;matrix6u;matrix7u;matrix4l;matrix5l;matrix6l];
matrix_raw=[matrix4u_raw;matrix5u_raw;matrix6u_raw;matrix7u_raw;matrix4l_raw;matrix5l_raw;matrix6l_raw];
origin=[4*ones(1,size(matrix4u,1)),5*ones(1,size(matrix5u,1)),6*ones(1,size(matrix6u,1)),7*ones(1,size(matrix7u,1)),...
    4*ones(1,size(matrix4l,1)),5*ones(1,size(matrix5l,1)),6*ones(1,size(matrix6l,1))];
%matrix7l appears to be subject to a mis-dissection. All cells originating
%in this sample will therefore be excluded from all future analysis

labels={'OB','SC','Li','LM','AL','PM','AM','RL','LP','LGN','Str','LeftVis','TeA/Ect'};


%% split spikes according to the mice they come from

%split into mice
xmouse4_tmp=x097(1:17);
xmouse5_tmp=x097(18:33);
xmouse6_tmp=x102(1:17);
xmouse7_tmp=x102(18:34);

%shuffle for consistent target areas, keeping ob and visual areas
% in order: ob li lm al pm am rl
m4=[4 6 9 10 11 12 13]; 
m5=[2 4 7 8 9 10 11 ];
m6=[3 5 8 9 10 11 12];
m7=[3 5 8 12 10 11 9];

xmouse(4,:)=xmouse4_tmp(m4);
xmouse(5,:)=xmouse5_tmp(m5);
xmouse(6,:)=xmouse6_tmp(m6);
xmouse(7,:)=xmouse7_tmp(m7);

xmouse=xmouse./repmat(xmouse(:,1),1,size(xmouse,2));

%% normalize data matrix using the spike in counts
load ../ZL102/threshold4/matrix_raw.mat
areas=[3,4,5,6,7,8]; %higher visual areas
matrix_raw=matrix_raw(:,[1,areas]); %remove everything but higher visual areas and OB negative control
matrix_raw=matrix_raw(max(matrix_raw(:,2:end),[],2)>10,:);

matrix_norm=[];
for m=4:7
    matrix_norm=[matrix_norm; matrix_raw(origin==m,:)./repmat(xmouse(m,:),sum(origin==m),1)];
end

save('matrix_norm.mat','matrix_norm','origin');



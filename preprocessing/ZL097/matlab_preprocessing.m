%% read in sequences and counts of data split into the different libraries

% %read in raw sequencing reads to determine BC+UMI+SSI read cutoff in
% shell preprocessing file
% 
% bcn=[97:132];
% for i=1:36
% data(i).rawcounts=dlmread(['ZL097BC',int2str(bcn(i)),'_counts.txt']);
% end
% save('rawdata.mat','data','-v7.3');
% 
% bcn=[97:132];
% figure;
% for i=1:36
%     loglog(data(i).rawcounts);
%     title(int2str(bcn(i)));
%     pause();
% end


bcn=[97:132];
for i=1:36
data(i).counts=dlmread(['ZL097',int2str(bcn(i)),'_counts.txt']);
data(i).reads=int8(char(textread(['ZL097',int2str(bcn(i)),'_seq.txt'],'%s')));
end

save('data.mat','data','-v7.3');


%% collapse using bowtie aligment results
positions1=[];

for i=1:36
positions1(i).x=dlmread(['bowtie',int2str(bcn(i)),'_2u_1.txt']);
positions1(i).y=dlmread(['bowtie',int2str(bcn(i)),'_2u_3.txt']);
clustermatrix1(i).C=sparse(positions1(i).x,positions1(i).y,1); %make a sparse matrix using the bowtie columns 1 and 3 as x and y coordinates for nonzero matrix entries
end

%find all connected components
graph=[];
for i=1:36
i
    [graph(i).S,graph(i).G]=graphconncomp(clustermatrix1(i).C,'Directed','false'); %find the connected graph components
end

%collapse all members of a connected component to its most abundant member,
%and so do error correction
for i=1:36
x=1:graph(i).S;
[tf,loc]=ismember(x,graph(i).G,'R2012a');
collapsedreads=data(i).reads(loc,:);
collapsedcounts=accumarray(graph(i).G',data(i).counts);%'
[corrected(i).counts2u,ix]=sort(collapsedcounts,'descend');
corrected(i).reads2u=collapsedreads(ix,:);
end


%remove reads containing homopolymers
minrunlength=7; % as 0.25^7*23=0.0014 or less than 1% of barcodes will have this by chance?
for i=1:36
    a=findhomopolymers(corrected(i).reads2u,minrunlength);
    corrected(i).freads=corrected(i).reads2u(~a,:);
    corrected(i).fcounts=corrected(i).counts2u(~a,:);
end


%% filter for match to virus library 

%load virus libary barcodes
load ../ZL077_viruslibary/collapsedslow

%check for overlap and keep only what matches the library
for i=1:36
d=[];    
d=ismember(corrected(i).freads,collapsed.reads,'rows');
data(i).BCseqff=corrected(i).freads(d,:);
data(i).BCcountsff=corrected(i).fcounts(d,:);
end
save('data.mat','data','-v7.3');
						   				   
		   
%% check out spike ins.

%first do the collapse.
						   
%read in non collapsed spike ins.
bcn=[97:132];

for i=1:36
    i
spikes(i).counts=dlmread(['ZL097spikes',int2str(bcn(i)),'_counts.txt']);
spikes(i).reads=int8(char(textread(['ZL097spikes',int2str(bcn(i)),'_seq.txt'],'%s')));
end
						   

%collapse using bowtie all-to-all alignment results
positions2=[];

for i=1:36
positions2(i).x=dlmread(['bowtiespikes',int2str(bcn(i)),'_2u_1.txt']);
positions2(i).y=dlmread(['bowtiespikes',int2str(bcn(i)),'_2u_3.txt']);
clustermatrix2(i).C=sparse(positions2(i).x,positions2(i).y,1); %make a sparse matrix using the bowtie columns 1 and 3 as x and y coordinates for nonzero matrix entries
end

graph2=[];
for i=1:36
i
    [graph2(i).S,graph2(i).G]=graphconncomp(clustermatrix2(i).C,'Directed','false'); %find the connected graph components
end

for i=1:36
x=1:graph2(i).S;
[tf,loc]=ismember(x,graph2(i).G,'R2012a');
collapsedreads=spikes(i).reads(loc,:);
collapsedcounts=accumarray(graph2(i).G',spikes(i).counts);%'
[spikes(i).counts2u,ix]=sort(collapsedcounts,'descend');
spikes(i).reads2u=collapsedreads(ix,:);
end


save('spikes.mat','spikes')



%% split into brains

%order data according to SSI choice
x=[1:36];

brain4=x(1:17);
brain5=x(18:33);

data4=[];
j=1;
for i=brain4
    data4(j).BCcountsff=data(i).BCcountsff;
    data4(j).BCseqff=data(i).BCseqff;
    j=j+1;
end

data5=[];
j=1;
for i=brain5
    data5(j).BCcountsff=data(i).BCcountsff;
    data5(j).BCseqff=data(i).BCseqff;
    j=j+1;
end



%% look at brain4 and build a barcode matrix

data=data4;

%define reference barcodes
refbarcodes_tmp=[];
for i=1:length(data)
refbarcodes_tmp=[refbarcodes_tmp;data(i).BCseqff];
end
refbarcodes=unique(refbarcodes_tmp,'rows');

barcodematrix=zeros(size(refbarcodes,1),size(data,2));%initiate the barcode matrix

j=1;
for i=1:size(data,2)
BCreads=data(i).BCseqff;
BCmolcounts=data(i).BCcountsff;
[ind,loc]=ismember(BCreads,refbarcodes,'rows');
barcodematrix(loc(loc~=0),i)=BCmolcounts(ind); 
end

matrix4=barcodematrix;
barcodes4=refbarcodes;
save('brain4.mat','matrix4','barcodes4')


%% look at brain 5 and build a barcode matrix

data=data5;

%define reference barcodes
refbarcodes_tmp=[];
for i=1:length(data)
refbarcodes_tmp=[refbarcodes_tmp;data(i).BCseqff];
end
refbarcodes=unique(refbarcodes_tmp,'rows');

barcodematrix=zeros(size(refbarcodes,1),size(data,2));%initiate the barcode matrix

j=1;
for i=1:size(data,2)
BCreads=data(i).BCseqff;
BCmolcounts=data(i).BCcountsff;
[ind,loc]=ismember(BCreads,refbarcodes,'rows');
barcodematrix(loc(loc~=0),i)=BCmolcounts(ind); 
end

matrix5=barcodematrix;
barcodes5=refbarcodes;
save('brain5.mat','matrix5','barcodes5')


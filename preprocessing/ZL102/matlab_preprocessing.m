function matlab_preprocessing
%% read in sequences and counts of data split into the different libraries
% Justus Kebschull

%% load data

% % check raw sequencing reads per BC+UMI+SSI combination; allows the
% % thresholding in the shell script
% bcn=[133:168];
% for i=1:36
% data(i).rawcounts=dlmread(['ZL102BC',int2str(bcn(i)),'_counts.txt']);
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

%load data
bcn=[133:168];
for i=1:36
	try
data(i).counts=dlmread(['ZL102',int2str(bcn(i)),'_counts.txt']);
data(i).reads=int8(char(textread(['ZL102',int2str(bcn(i)),'_seq.txt'],'%s')));
catch
	end
end

save('data.mat','data','-v7.3');


%% error correct the raw sequencing reads

%collapse using bowtie 
positions1=[];
graph=[];
for i=1:36
	try
positions1(i).x=dlmread(['bowtie',int2str(bcn(i)),'_2u_1.txt']);
positions1(i).y=dlmread(['bowtie',int2str(bcn(i)),'_2u_3.txt']);
clustermatrix1(i).C=sparse(positions1(i).x,positions1(i).y,1); %make a sparse matrix using the bowtie columns 1 and 3 as x and y coordinates for nonzero matrix entries

%find all connected components
[graph(i).S,graph(i).G]=graphconncomp(clustermatrix1(i).C,'Directed','false'); %find the connected graph components
x=1:graph(i).S;
%collapse all barcodes within a connected component to its most abundant
%member
[tf,loc]=ismember(x,graph(i).G,'R2012a');
collapsedreads=data(i).reads(loc,:);
collapsedcounts=accumarray(graph(i).G',data(i).counts);%'
[corrected(i).counts2u,ix]=sort(collapsedcounts,'descend');
corrected(i).reads2u=collapsedreads(ix,:);
catch
	end
end


%remove reads containing homopolymers
minrunlength=7; % as 0.25^7*23=0.0014 or less than 1% of barcodes will have this by chance?
for i=1:36
		try
    a=findhomopolymers(corrected(i).reads2u,minrunlength);
    corrected(i).freads=corrected(i).reads2u(~a,:);
    corrected(i).fcounts=corrected(i).counts2u(~a,:);
	catch 
		end
end


%% filter for match to virus library 

%load virus libary barcodes
load ../ZL077_viruslibrary/collapsedslow

%check for overlap and keep only what matches the library
for i=1:36
			try
d=[];    
d=ismember(corrected(i).freads,collapsed.reads,'rows');
data(i).BCseqff=corrected(i).freads(d,:);
data(i).BCcountsff=corrected(i).fcounts(d,:);
catch
	end
		
end
save('data.mat','data','-v7.3');
						   				   

				   
%% check out spike ins.

%first do the collapse.
						   
%read in non collapsed spike ins.
bcn=[133:168];

for i=1:36
	try
    i
spikes(i).counts=dlmread(['ZL102spikes',int2str(bcn(i)),'_counts.txt']);
spikes(i).reads=int8(char(textread(['ZL102spikes',int2str(bcn(i)),'_seq.txt'],'%s')));
catch
	end
end
						   
						   


%collapse using bowtie alignments as above
positions2=[];
graph2=[];

for i=1:36
	try
positions2(i).x=dlmread(['bowtiespikes',int2str(bcn(i)),'_2u_1.txt']);
positions2(i).y=dlmread(['bowtiespikes',int2str(bcn(i)),'_2u_3.txt']);
clustermatrix2(i).C=sparse(positions2(i).x,positions2(i).y,1); %make a sparse matrix using the bowtie columns 1 and 3 as x and y coordinates for nonzero matrix entries

i
    [graph2(i).S,graph2(i).G]=graphconncomp(clustermatrix2(i).C,'Directed','false'); %find the connected graph components

x=1:graph2(i).S;
[tf,loc]=ismember(x,graph2(i).G,'R2012a');
collapsedreads=spikes(i).reads(loc,:);
collapsedcounts=accumarray(graph2(i).G',spikes(i).counts);%'
[spikes(i).counts2u,ix]=sort(collapsedcounts,'descend');
spikes(i).reads2u=collapsedreads(ix,:);
catch
	end
end


save('spikes.mat','spikes')


%% split data into brains

%order data according to SSI choice
x=[1:36];

brain6=x(1:17);
brain7=x(18:34);

data6=[];
j=1;
for i=brain6
    data6(j).BCcountsff=data(i).BCcountsff;
    data6(j).BCseqff=data(i).BCseqff;
    j=j+1;
end

data7=[];
j=1;
for i=brain7
    data7(j).BCcountsff=data(i).BCcountsff;
    data7(j).BCseqff=data(i).BCseqff;
    j=j+1;
end



%% look at brain6 and build a barcode matrix
data=data6;


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

matrix6=barcodematrix;
barcodes6=refbarcodes;
save('brain6.mat','matrix6','barcodes6')


%% look at brain7 and build a barcode matrix
data=data7;


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

matrix7=barcodematrix;
barcodes7=refbarcodes;
save('brain7.mat','matrix7','barcodes7')

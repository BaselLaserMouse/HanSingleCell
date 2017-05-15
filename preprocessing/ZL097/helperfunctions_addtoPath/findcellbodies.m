function [cellbodylocation,data_filtered,barcodes_filtered,keep]=findcellbodies(data,ratiocutoff,barcodes)


%find max
[firstmax,loc]=max(data,[],2);
%find second max
data_max_tmp=data;
for i=1:length(data);
data_max_tmp(i,loc(i))=0;
end
[secondmax]=max(data_max_tmp,[],2);

maxratio=secondmax./firstmax; %take ratio
% figure; plot(firstmax,secondmax,'.');
% figure;hist(maxratio,1000);

%keep only those barcodes for which ratio is smaller than specified amount

data_filtered=data(maxratio<ratiocutoff,:);
barcodes_filtered=barcodes(maxratio<ratiocutoff,:);

keep=maxratio<ratiocutoff;

%locate position of cell bodies
cellbodylocation=zeros(size(data_filtered));
[m,cell_loc]=max(data_filtered,[],2);
for i=1:length(data_filtered)
    cellbodylocation(i,cell_loc(i))=1;
end


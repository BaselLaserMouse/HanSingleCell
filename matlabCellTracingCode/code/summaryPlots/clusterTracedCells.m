function varargout=clusterTracedCells(pltData,excludeAreas,linkageThresh,numNodes)
% hierarchical clustering of traced neurons by brain area
%
% function stats=clusterTracedCells(pltData,excludeAreas,linkageThresh)
%
% Purpose
% hierarchical clustering of traced neurons by brain area
%
% Inputs
% pltData - a structure that is the output of pointsByAreaPlot.
% excludeAreas - a cell array of strings indicating which areas to
%               exclude from the analysis. These are the long 
%				names. The string 'Primary visual area' will match
%				'Primary visual area, layer 4', 'Primary visual area, layer 6', etc...
%				NOTE: this just does a string-based match. No reference
% 				to the hierarchical organisation of the tree is made. 
% linkageThresh - The threshold at which to highligh different parts of the plotted 
%                 dendrogram (absent by default).
% numNodes		- How many terminal nodes to put in the dendrogram. By default, each 
% 				  cell is one node
%
%
% Rob Campbell - 2016


if nargin<2
	excludeAreas=[];
end

if nargin<3 | isempty(linkageThresh)
	linkageThresh='default';
end

if nargin<4 | isempty(numNodes)
	numNodes=0; %if zero the number of nodes is not restricted
end


%Exclude if needed
[pltData,excluded] = filterPlotAreasFromPltData(pltData,excludeAreas);

d=pltData.dataMat';

dist = pdist(d,'cosine');
clustTreeEuc = linkage(dist,'single');


[h,nodes,perm] = dendrogram(clustTreeEuc,numNodes,'ColorThreshold',linkageThresh);
set(h,'LineWidth',2)
if numNodes==0
	set(gca,'XTick',1:size(pltData.dataMat,2), 'XTickLabel',pltData.cellIDs(perm),'XTickLabelRotation',45)
	xlabel('Sample')
else
	xlabel('Group')
end

ylabel('Linkage')


if nargout>0
	stats.clustTreeEuc=clustTreeEuc;
	stats.perm=perm;
	stats.dist=dist;
	stats.nodes=nodes;
	stats.excluded=excluded;
	stats.notes='';
	varargout{1} = stats;
end


function neuriteTrees=getTreeInArea(data,areaID,doDiagnosticPlot)
% Return the sub-trees of a neurite tree that have terminal nodes located within ARA area having index areaID
%
%
%	function getTreeInArea(data,areaID)
%
% data is the output of getLeavesInARA


if nargin<3
	doDiagnosticPlot=false;
end



if doDiagnosticPlot
	clf
end
labels=getAllenStructureList;




L = data.neuriteTree.findleaves;
for ii=2:length(L) %The first one is the root node ID
	lastNodeIndex(ii) = lastNodeInArea(data,L(ii),doDiagnosticPlot,labels);

end



	function lastNodeIndex = lastNodeInArea(data,nodeID,doDiagnosticPlot,labels)
		nTree = data.neuriteTree;
		lastNodeAreaIndex = data.allNodeInd(nodeID);

		nodeIndToRoot=nTree.pathtoroot(nodeID);
		lastNodeIndex=nan;

		if doDiagnosticPlot
			inds = data.allNodeInd(nodeIndToRoot);
			%Skip if all inds are the same (then it's an axon in the same area as the soma)
			if length(unique(inds))==1
				return
			end
			plot(inds,'ok-','MarkerFaceColor',[1,1,1]*0.75)
			
			indLabels =unique(inds);
			for kk=1:length(indLabels)
				indLabelsText{kk} = sprintf('%s (%d)',labels.name{labels.id == indLabels(kk)},indLabels(kk)) ;
			end
			set(gca,'YTick',indLabels,'YTickLabel',indLabelsText)
			grid on
			title( sprintf('path %d/%d',ii,length(L)) )
			drawnow
			pause
		end


	end


end
function axonTree2Dplots(neuriteTree,highlightNode,allOfSubTree)
% Make sagittal, coronal, and transverse plots of the neurite tree
%
% function axonTree2Dplots(neuriteTree,highlightNode,allOfSubTree)
%
% Purpose
% This is mainly for diagnostic purposes. Does not overlay brain. 
% Just shows the tree shape in three subplots.
%
% Inputs
% neuriteTree - a tree object produce by MaSIV.
% highlightNode - an optional index which, if present, highlights this node on the plot
% allOfSubTree - [0 by default] if 1 all nodes in the sub-tree after highlightNode are marked
%
%
% Rob Campbell - Basel 2016



if isempty(neuriteTree)
	fprintf('neuriteTree is empty. aborting\n')
	return
end

if nargin<2
	highlightNode=[];
end

if nargin<3
	allOfSubTree=0;
end


clf

subplot(2,2,1)
plotTree([1,2])
return
subplot(2,2,2)
plotTree([2,3])

subplot(2,2,3)
plotTree([1,3])





function plotTree(dims)
	%Get the data from the tree
	segments = neuriteTree.getsegments;

	hold on
	for ii=1:length(segments)
		theseNodes = neuriteTree.Node(segments{ii});
		xyz = ones(length(theseNodes),3);

		for jj=1:length(theseNodes)  %Because the MaSIV trees are different to the re-imported trees made by exportedCSV2tree
			if isa(theseNodes{1},'neuriteTracerNode')                   
				xyz(jj,:) = [theseNodes{jj}.xVoxel,theseNodes{jj}.yVoxel,theseNodes{jj}.zVoxel];
			else
				xyz(jj,:) = theseNodes{jj};
			end
		end

		plot(xyz(:,dims(1)), xyz(:,dims(2)), '-', 'linewidth', 1, 'color', 'r');
	end %for ii=1:length(segments)

	hN=highlightNode;

	if ~isempty(hN) & hN<=length(neuriteTree.Node)
		if isa(neuriteTree.Node{1},'neuriteTracerNode') 
			H = [neuriteTree.Node{hN}.xVoxel, neuriteTree.Node{hN}.yVoxel, neuriteTree.Node{hN}.zVoxel];
		else
			H = [neuriteTree.Node{hN}];
		end
		if ~isempty(H)
			plot(H(dims(1)), H(dims(2)), '*b')
		else
			fprintf('No node %d found\n',highlightNode)
		end
		ST = neuriteTree.subtree(hN);

		if allOfSubTree
			xyz = ones(length(ST.Node),3);
			for ii=1:length(ST.Node)
			if isa(neuriteTree.Node{1},'neuriteTracerNode') 
				xyz(ii,:) = [ST.Node{ii}.xVoxel,ST.Node{ii}.yVoxel,ST.Node{ii}.zVoxel];
			else
				xyz(ii,:) = ST.Node{ii};
			end
				plot(xyz(ii,dims(1)), xyz(ii,dims(2)), '.k')
			end
		end %if allOfSubTree

	end %if ~isempty(hN) & hN<=length(neuriteTree.Node)

	hold off

	axis square
	box on
	grid on
	drawnow
	end %close plotTree

end %close main function body
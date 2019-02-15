function varargout = compareTraces(traces)
% compare neurite traces
%
%
%
% Purpose
% Return a simple comparison of number of neurite trees 
%
% Inputs
% traces - a cell array of neurite trees to compare. each cell should contain 
%		a different trace of the same cell.
%
%
%
%
% Rob Campbell - Basel 2016


if ~iscell(traces)
	error('traces should be a cell array');
end

for ii=1:length(traces)
	fprintf('%d: %d leaves, %d segments\n',ii,length(traces{ii}.findleaves),length(traces{ii}.getsegments))
end


clf

subplot(2,2,1)
plotTrees(traces,'x','y')

subplot(2,2,2)
plotTrees(traces,'x','z')

subplot(2,2,3)
plotTrees(traces,'y','z')


function plotTrees(traces,dimA,dimB)

	hold on
	cols=parula(length(traces));
	for ii=1:length(traces)
		t=traces{ii};

		segments=t.getsegments;
		for jj=1:length(segments)
			theseNodes = t.Node(segments{jj});
			data = segment2points(theseNodes,dimA,dimB);
			plot(data(:,1), data(:,2), '-', 'color', cols(ii,:))
		end %jj

	end

	%soma
	plot(t.Node{1}.([dimA,'Voxel']),t.Node{1}.([dimB,'Voxel']),'or')

	grid on 
	box on
	set(gca,'Color',[1,1,1]*0.8)
	drawnow


function data = segment2points(thisSegment,dimA,dimB)
	data = ones(length(thisSegment),2);

	for ii=1:length(thisSegment)
		data(ii,1) = thisSegment{ii}.([dimA,'Voxel']);
		data(ii,2) = thisSegment{ii}.([dimB,'Voxel']);
	end
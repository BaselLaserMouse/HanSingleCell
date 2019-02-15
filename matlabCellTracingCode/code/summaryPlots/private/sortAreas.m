function out = sortAreas(pltData,sortByStrength)
% pointsByAreaPlot helper function - sorts Allen Reference Atlas areas for plotting
%
% function out = sortAreas(pltData,sortByStrength)
%
% Purpose
% Nicely sort brain area names for area plot (pointsByAreaPlot)
% sorts the pltData variable produced by pointsByAreaPlot and called 
% within this function. By number of connections and groups areas by layer
% 
%
% Inputs
% pltData - made by pointsByAreaPlot
% sortByStrength - [true by default] sorts areas by connection strength.
%
% Outputs
% out - a sorted version of pltData
%
%
% Rob Campbell - Basel 2015



if nargin<2
	sortByStrength=true;
end

%Find all names that indicate a name refers to a particular layer within a brain (probably cortical) region 
ind=regexp(pltData.areaNamesInSamples,'[lL]ayer \d');

layers = {};
nonLayerInd=[]; %Keeps track of which brain regions don't have a layer name associated with them
for ii=1:length(ind)
	indInString = ind{ii};
	
	if isempty(indInString)
		nonLayerInd(length(nonLayerInd)+1)=ii;
		continue
	end
	layers = [layers,pltData.areaNamesInSamples{ii}(1:indInString-1)];

end

uniqueLayerNames = unique(layers)';

if length(layers)+length(nonLayerInd) ~= length(pltData.areaNamesInSamples)
	error('Appear to be missing areas (areas with layers plus areas without layers do not sum to original number of areas)')
end


%Loop through the area names and sort into layers where possible
groupedAreas={};
for ii=1:length(uniqueLayerNames)

	%must match non-unique entries only. This ensures that the parent structure (without layers attached to the name)
	%does not appear more than once. This is a horrible hack and we need to get this shit working nicely using the 
	%tree hiearchy. 
	f=strmatch(uniqueLayerNames{ii},pltData.areaNamesInSamples); 
	e=strmatch(uniqueLayerNames{ii},pltData.areaNamesInSamples,'exact');	
	if ~isempty(e)
		f(f==e)=[]; %remove the extra (non-layer) entry from the group
	end

	[~,ind]=sort(pltData.areaNamesInSamples(f));

	f=f(ind); %The sorted indexes from one area, which we will treat as one entity
	groupedAreas{length(groupedAreas)+1}=f;

end


%we want to sort the grouped areas (those with sub-layers) by the number of hits
groupedAreas=[groupedAreas,num2cell(nonLayerInd)];

T=cellfun(@(x) sum(pltData.allAreas(x,2)), groupedAreas);

if sortByStrength
	[~,ind] = sort(T,'descend');
	groupedAreas = groupedAreas(ind);
end

out.allAreas = [];
out.dataMat = [];
out.areaNamesInSamples={};
out.groupingPoints=[]; %The breakpoints for where groups of areas end and start



for ii=1:length(groupedAreas)
	ind = groupedAreas{ii};

	out.allAreas=[out.allAreas; pltData.allAreas(ind,:)];
	out.dataMat=[out.dataMat; pltData.dataMat(ind,:)];
	out.areaNamesInSamples = [out.areaNamesInSamples; pltData.areaNamesInSamples(ind)];

	if length(ind)>1
		n=length(out.areaNamesInSamples);
		out.groupingPoints(length(out.groupingPoints)+1,1:2)=[n-length(ind),n];
	end
end


out.groupingPoints = unique(out.groupingPoints(:));
out.groupingPoints(out.groupingPoints==0)=[];


%Report if the plot matrix has changed size
if ~isequal(size(pltData.dataMat),size(out.dataMat))
	fprintf('\n ** %s: dataMat has changed size! from %dx%d to %dx%d **\n\n', mfilename, size(pltData.dataMat),size(out.dataMat))

	%If the first dimension is now smaller than figure out which area names are missing
	if size(pltData.dataMat,1) > size(out.dataMat,1)
		fprintf('Now missing the following brain areas:\n')
		for ii=1:length(pltData.areaNamesInSamples)
			if isempty(strmatch(pltData.areaNamesInSamples{ii}, out.areaNamesInSamples,'exact'))
				disp(pltData.areaNamesInSamples{ii})
			end
		end
	end

end


%Add missing fields
out.dataStructure = pltData.dataStructure;
out.cellIDs = pltData.cellIDs;
out.somaLocation=pltData.somaLocation;
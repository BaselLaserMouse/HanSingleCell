function varargout=overlayCellsOnThreeProjections(cellFnames,projectedAreas,displayProps)
% Add one more cells to all three projection directions as a single figure
%
% function [H,areas]=overlayCellsOnThreeProjections(cellFnames,projectedAreas,displayProps)
%
% Purpose
% Make a publication-ready plot of one or more cells overlaid onto the three projections
% Plots include midline, cortical boundary, and V1. These can be altered by modifying the
% structure of function handles returned by this function. 
%
% NOTE: this function assumes we're working with the 25 micron ARA and it requires you
%       to be in the project root directory. e.g. the hanyu folder on the server. You can 
%       also run the function from a different directory, so long as you define the ovDir
%       argument using the absolute path to the outlines and also the absolute path to the 
%       cellFnames.
%
%
% Inputs
% cellFnames - A cell array of strings listing relative or absolute paths
%             to traced cell file names. If empty, no cells are plotted.
% projectedAreas - a) The output of generateProjections. 
%                  b) The output of this function. If so, any structure in the third argument
%                     is ignored, all existing neurons are removed and new ones are added.
%                     This should speed up drawing of complicated boundaries, and make user 
%                     modification of existing boundaries more flexible.
% displayProps - [optional] A cell array that is used to control the appearance of the plot. 
%                displayProps defines the plotting properties of each neuron:
%                - If it contains one cell, it applies this to all neurons. (use {'color','g'} not {'g'})
%                - If it's a cell array of cell arrays then it must have length(cellFnames) 
%                  and it sets the properties of each neuron in turn.
%                - If a structure, it must have the format for the second output of
%                  for instance, brainAreaNames.visualAreas. This will colour-code
%                  the listed areas and also the cell itself. 
%
%
% Outputs
% H - A structure of plot handles broken down by brain orientation
% areas - A structure of plot handles grouped by outline name, so all outlines in all 
%         plots can be changed at the same time. TODO
%
%
%
% Examples:
% In all cases run from experiment root directory and all examples need the following
% to have been executed before running them:
% L=listAllTracedCells;
% load ../ReferenceAtlas/ARA_CCFv3/ARA_25_micron_mhd/projections/brainAreaProjections.mat %or make the projections
%
%
%
% One - plotting one cell
%  overlayCellsOnThreeProjections(L('YH121_01'),projections) %one cell in red
%  overlayCellsOnThreeProjections(L('YH121_01'),projections,{'color','g'}) %or make it green at plot time
%
%
% Two - plotting many cells in different colours
%  cellPaths =  {L('YH265_01'),L('YH183_01'),L('YH290_01')};
%  props = {{'color','g'}, {'color','b'}, {'color','r','linewidth',2}};
%  overlayCellsOnThreeProjections(cellPaths,projections,props) %make cells different colors
%
%
% Three - color-coding the visual areas
%  [n,c]=brainAreaNames.visualAreas;
%  cellPaths =  {L('YH265_01'),L('YH183_01'),L('YH290_01')};
%  overlayCellsOnThreeProjections(cellPaths,projections,c)
%
%
% Four - preserving formatting 
% clf
% [b,c]=brainAreaNames.visualAreas; 
% H=overlayCellsOnThreeProjections(L('YH265_01'),projections,c);
% set(H.transverse.areas(1:3:end),'edgecolor','r')
% H=overlayCellsOnThreeProjections(L('YH121_01'),H,c);  %re-run with a different cell
%
%
%
% Rob Campbell - Basel 2016
%
%
% See also:
% overlayCellOnARA
% listAllTracedCells
% overlayManyCellsOnOutline



% Process input args
if nargin<1 
    help(mfilename)
    return
end

if nargin<2
    error('You need to supply both the cell names and the projectedAreas structure.')
end

if ~isstruct(projectedAreas)
    error('projectedAreas must be a structure produced by generateProjections or the output of this function')
end

if isfield(projectedAreas,'axesSagittal')
    %Then we remove already existing neurons and overlay new neurons onto the existing area boundaries
    H=projectedAreas;
    delete(H.coronal.neuriteTree{:}), H.coronal.neuriteTree=[];
    delete(H.transverse.neuriteTree{:}), H.transverse.neuriteTree=[];
    delete(H.sagittal.neuriteTree{:}), H.sagittal.neuriteTree=[];

    delete(H.coronal.somaMarker{:}), H.coronal.somaMarker=[];
    delete(H.transverse.somaMarker{:}), H.transverse.somaMarker=[];
    delete(H.sagittal.somaMarker{:}), H.sagittal.somaMarker=[];

    delete(H.coronal.highlightedLeaves{:}), H.coronal.highlightedLeaves=[];
    delete(H.transverse.highlightedLeaves{:}), H.transverse.highlightedLeaves=[];
    delete(H.sagittal.highlightedLeaves{:}), H.sagittal.highlightedLeaves=[];
    useExisting=true;
else
    useExisting=false;
end

if nargin<3
    displayProps=[];
end


if ischar(cellFnames)
    cellFnames = {cellFnames};
end

if isempty(cellFnames)
    cellFnames={};
end
if ~iscell(cellFnames)
    error('cellFnames is not a cell array')
end


cellColors=[];
areaHighlights=[];
if ~isempty(displayProps)
    if iscell(displayProps)
        if length(displayProps) ~= length(cellFnames) && iscell(displayProps{1})
            fprintf('You asked for %d cells but %d cell props. The numbers should be equal to color-code cells differently.\n',...
                length(cellFnames), length(displayProps))
            displayProps=[];
        end
        cellColors=displayProps;
    end

    if isstruct(displayProps)
        areaHighlights=displayProps;
        %Then we will also need the pointsInARA data for this neuron 

        for ii=1:length(cellFnames)
            thisFname = strrep(cellFnames{ii},'.csv','_pointsInARA.mat');
            if ~exist(thisFname)
                fprintf('%s can not find file %s. Skipping neurite color-coding\n',mfilename,thisFname)
                pointsInARA(ii)=struct;
                continue
            end

            load(thisFname)

            %This will do as a first pass, we may later want to load the xylem structure for more flexibility. 
            areaHighlights.pointsInARA(ii) = ARApoints.pointsInARA.rawSparseData; 
        end %for ii=1:length(cellFnames)

    end

end


% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
if ~useExisting
    clf 
    set(gcf,'PaperPosition',[0,0,12,12]) %Set the figure size to a fixed value to make saving more consistent
end

for ii=1:length(cellFnames)
    cellFnames{ii} = strrep(cellFnames{ii},'.csv','_pointsInARA.mat');
end

%Make the transverse plot
if ~useExisting
    H.axesTransverse=axes('position',[0.01,0.02,0.47,1],'tag','transverse');
else
    axes(H.axesTransverse)
end
H.transverse = makePlot(cellFnames,projectedAreas,areaHighlights,'transverse',1);
axis xy

%Make the coronal plot
if ~useExisting
    H.axesCoronal=axes('position',[0.52,0.46,0.47,0.47],'Tag','coronal');
else
    axes(H.axesCoronal)
end
H.coronal = makePlot(cellFnames,projectedAreas,areaHighlights,'coronal',0);

%Make the sagittal plot
if ~useExisting
    H.axesSagittal=axes('position',[0.52,0.1,0.47,0.47],'tag','sagittal');
else
    axes(H.axesSagittal)
end
H.sagittal = makePlot(cellFnames,projectedAreas,areaHighlights,'sagittal',0);


%Set the figure size to a default value to make saving more consistent
set(gcf,'PaperPosition',[0,0,30,29],'Name','Cell projections') 





%Format some of the areas
formatPlot(H.transverse,cellColors)
formatPlot(H.coronal,cellColors)
formatPlot(H.sagittal,cellColors)

%Optionally return the handles as a structure
if nargout>0
    varargout{1}=H;
end

if nargout>1
    f=unique([fields(H.coronal); fields(H.transverse); fields(H.sagittal)]);

    orientations = fields(H);

    areas = struct;

    for thisField=1:length(f)
        if strcmp(f{thisField},'neuronGroup')
            continue
        end

        areas.(f{thisField}) = [];

        for ii = 1:length(orientations)
            if isfield(H.(orientations{ii}),(f{thisField}))
                areas.(f{thisField}) = [areas.(f{thisField}), H.(orientations{ii}).(f{thisField})];
            end
        end
    end

    varargout{2} =areas;
end






% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
function H=makePlot(cellFnames,projectedAreas,areaHighlights,orientation,swapAxes)
    %Default  line properties
    outlineProps = {'-','linewidth',2,'color',[0.5,0.5,0.5]};
    cellProps = {'-','linewidth',1,'color', 'k'};
    sectionNames={'transverse','sagittal','coronal'};


    set(gcf,'PaperPosition',[0,0,12,12]) %Set the figure size to a fixed value to make saving more consistent
    axis off ij equal

    hold on 
        
    %add area outlines to figure
    switch orientation
        case 'transverse'
            dimsToPlot = [1,2];
        case 'sagittal'
            dimsToPlot = [2,1];
    case 'coronal'
            dimsToPlot = [2,1];
    end

    if swapAxes
        dimsToPlot = fliplr(dimsToPlot);
    end

    %PLOT PROJECTIONS
    ind=strmatch(orientation,sectionNames);

    if ~isfield(projectedAreas,orientation) %Then we don't modify anything
        H=projectionPlotter(projectedAreas(ind),areaHighlights);
    else
        H=projectedAreas.(orientation);
    end

    %add cells to figure
    switch orientation
        case 'transverse'
            dimsToPlot = [1,2];
        case 'sagittal'
            dimsToPlot = [1,3];
    case 'coronal'
            dimsToPlot = [2,3];
    end

    if swapAxes
        dimsToPlot = fliplr(dimsToPlot);
    end

    for ii=1:length(cellFnames)
        thisItem = cellFnames{ii};
        if isstr(thisItem)
            [H.neuriteTree{ii},H.somaMarker{ii},H.highlightedLeaves{ii}]=plotThisTree(thisItem,dimsToPlot,cellProps,areaHighlights);
        elseif iscell(thisItem)
            for kk=1:length(thisItem)
                [H.neuriteTree{ii}{kk},H.somaMarker{ii}{kk},H.highlightedLeaves{ii}]=plotThisTree(thisItem{kk},dimsToPlot,cellProps,areaHighlights);
            end
        else
            error('Cell name is neither a cell array or a string. quitting.')
        end
    end

    hold off



%---------------------------------------------------------------------------
function [neuriteTree,somaMarker,highlightedLeaves]=plotThisTree(dataFname,dimsToPlot,lineProps,areaHighlights)
    % Loads a neurite tree from a CSV file then plots it
    if ~exist(dataFname,'file')
        fprintf('Skipping file %s. Can not find it\n',dataFname)
        neuriteTree=[];
        somaMarker=[];
        highlightedLeaves=[];
        return
    end

    load(dataFname)
    data = exportedCSV2tree(ARApoints.pointsFname);
    plotTree=1;

    if plotTree

        segments = data.getsegments;
        for ii=1:length(segments)
            theseNodes = data.Node(segments{ii});

            theseData = ones(length(theseNodes),3);
            for jj=1:length(theseNodes)
                theseData(jj,:) = theseNodes{jj};
            end

            neuriteTree(ii)=plot(theseData(:,dimsToPlot(1)), theseData(:,dimsToPlot(2)), lineProps{:});
        end
    else
        U=upsampleNeuriteTree(data,25);
        for ii=1:10:length(U)
            neuriteTree(1)=plot(U(ii,dimsToPlot(1)), U(ii,dimsToPlot(2)), '.r');
        end
    end

    %Plot the soma
    soma = data.Node{1};
    somaMarker=plot(soma(:,dimsToPlot(1)), soma(:,dimsToPlot(2)), 'ow','MarkerFaceColor','m','MarkerSize',8);

    %Labels:
    % bright termination = square
    % fading termination = triangle
    % grey matter = blue
    % white matter = green

    %Overlay points highlighted by the person who did the tracing
    oTrace=ARApoints.origTrace;
    highlightedLeaves=[];
    wInds=aratools.utils.whiteMatterInds;

    for ii=1:length(oTrace.Node)
        N=oTrace.Node{ii};

        if ~isempty(N.data)
            if isfield(N.data,'nodeType')
                nType=N.data.nodeType;

                pltSettings={};
                if regexp(nType,'bright','once')
                    pltSettings=[pltSettings, {'Marker','s'}];
                elseif regexp(nType,'fading','once')
                    pltSettings=[pltSettings, {'Marker','^'}];
                end

                if isempty(pltSettings)
                    %There are small mis-matches between the downsampled traces and the ones in X. e.g. an added node on a final
                    %branch because someone added a label. So we need this here.
                    continue
                end

                if any(wInds==ARApoints.pointsInARA.rawSparseData.ARAindex(ii)) %Then it's white matter
                    pltSettings=[pltSettings, {'MarkerFaceColor',[0.33,0.9,0.33]}];
                else
                    pltSettings=[pltSettings, {'MarkerFaceColor',[0.5,0.5,1]}];
                end
                coords=data.Node{ii};
                highlightedLeaves(end+1)=plot(coords(:,dimsToPlot(1)), coords(:,dimsToPlot(2)), 'w', 'MarkerSize',8, pltSettings{:});

            end
        end
    end




function H=projectionPlotter(projectionStructure,highlightAreas)

    % Plots the area boundaries
    H.areas=[];
    n=height(projectionStructure.structureList);

    for ii = 1:n

        B = projectionStructure.structureList.areaBoundaries{ii}; %Collect the border data for this area

        thisName=projectionStructure.structureList.name{ii};

        if ~isempty(highlightAreas) && ~isempty(strmatch(thisName,highlightAreas.areaNames))
            %We highlight the area
            thisName=strrep(thisName,' ','_');
            for k = 1:length(B)
                thisBoundary = B{k};
                %skip really tiny boundaries. They are artifacts
                if length(thisBoundary)<15
                    continue
                end
                tmp=patch(thisBoundary(:,2), thisBoundary(:,1),highlightAreas.areas(thisName));
                set(tmp, 'EdgeColor', highlightAreas.areas(thisName), ...
                    'FaceColor', highlightAreas.areas(thisName), ...
                    'LineWidth',1, ...
                    'Tag', strrep(thisName,'_',' ') );

                H.areas=[H.areas,tmp];
            end
        else
            %Just draw the regular boundary
            for k = 1:length(B)
                thisBoundary = B{k};
                %skip really tiny boundaries. They are artifacts
                if length(thisBoundary)<10
                    continue
                end
                tmp = patch(thisBoundary(:,2), thisBoundary(:,1), [0.98,0.97,0.97], 'EdgeColor', [1,1,1]*0.75, 'linewidth',0.5);
                set(tmp, 'Tag', 'AREA BOUNDARY');
                H.areas=[H.areas, tmp];
            end
        end
    end



function formatPlot(H,displayProps)
    if isfield(H,'primary_visual')
        %set(H.primary_visual,'LineStyle','--') 
    end

    if isempty(displayProps)
        return 
    end

    if ~iscell(displayProps{1}) %apply this one set of props to all
        set([H.neuriteTree{:}],displayProps{:})
        return
    end

    for ii=1:length(displayProps)
        set(H.neuriteTree{ii},displayProps{ii}{:})
    end



%% Figure to plot heatmaps of MAPseq projection strengths onto brain space
function [H]=make_nice_MAPseqplots3(atlas,loadings,marker,width)
if nargin < 3
    marker=3;
    width=1;
end


%load projections/brainAreaProjections.mat

%set areas
areas={'Laterointermediate area','Lateral visual area','Anterolateral visual area','posteromedial visual area',...
    'Anteromedial visual area','Rostrolateral visual area'};

%load colormap
load heatmap_colormap.mat
area_colors_template=cmap2;
area_colors=area_colors_template(ceil(loadings.*63)+1,:);


%plot figure;
H=aratools.projectAtlas.highlightAreaPlot(atlas(1),areas);
xlim([275 400])
ylim([110 243])


for i=1:length(areas)
    try
set(H.hLight(strmatch(areas(i),H.plottedNames),:),'FaceColor',area_colors(i,:));
    catch
    end
end

% add a cross to mark the injection spots
hold on
plot(350,150,'x','MarkerSize',marker,'LineWidth',width,'Color',area_colors_template(end,:))
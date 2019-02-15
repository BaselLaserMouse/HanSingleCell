function dataStruct=buildThumbnails(X)
% Build hemifield (thumbnail) images of each cell in the database. 
%
% function dataStruct=scc.diagnostic.buildThumbnails(X)
%
%
% Purpose
% Call from project root dir.
%
% Inputs
% X is the xylem data object containing all the traced cells 


if nargin<1
    help(['scc.diagnostic.',mfilename])
    return
end

tic

load('../ReferenceAtlas/ARA_CCFv3/ARA_25_micron_mhd/projections/brainAreaProjections.mat')


if ~exist('DiagnosticPlots','dir')
    error('Expected to find DiagnosticPlots directory. It is not there')
end

targetDir='./DiagnosticPlots/thumbnailsDatabase';
imageDir=fullfile(targetDir,'images');
if ~exist(imageDir,'dir')
    fprintf('Making directory %s\n', imageDir);
    mkdir(imageDir)
end




[visualAreaNames,colorMap]=brainAreaNames.visualAreas;


clf
set(gcf,'Color','w','Renderer','Painters')


% Loop through and make each cell plot
H=overlayCellsOnThreeProjections(X.data(1).pointsFname,projections,colorMap);

%Overlay sub-cortical areas present on this neuron
ii=1;
S=overlaySelectAreas(H,X.data(ii).pointsInARA.leaves.ARAindex);
thisFname = fullfile(imageDir,X.data(ii).details.cellID);
makeAndSaveThumbnail(H,thisFname)

for ii=1:length(X.data)

    fprintf('%d/%d. %s\n', ii, length(X.data), X.data(ii).pointsFname)
    thisFname = fullfile(imageDir,X.data(ii).details.cellID);
    delete(S)
    H=overlayCellsOnThreeProjections(X.data(ii).pointsFname,H,colorMap);
    S=overlaySelectAreas(H,X.data(ii).pointsInARA.leaves.ARAindex);
    T=text(350,380, strrep(X.data(ii).details.cellID,'_','\_'), 'parent',H.axesTransverse,...
        'FontSize',18,'FontWeight','Bold');

    makeAndSaveThumbnail(H,thisFname)

    delete(T)
end




    function makeAndSaveThumbnail(H,thisFname)

        thumbFig=figure;
        ax=copyobj(H.axesTransverse,thumbFig);
        set(ax,'Position',[0,0,1,1])
        xlim([230,435])
        ylim([105,390])

        %Add a scale bar
        hold on
        plot([300,340], [380,380], '-k', 'LineWidth', 4)
        hold off
        
        fprintf('Saved image to %s... ',thisFname)
        thisFname=[thisFname,'.eps'];
        print(gcf,'-depsc',thisFname)
        if exist('epsclean','file')
            epsclean(thisFname)
        end
        fprintf('DONE.\n');

        delete(thumbFig)


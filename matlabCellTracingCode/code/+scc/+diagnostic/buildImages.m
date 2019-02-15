function dataStruct=buildImages(X)
% Build nice images of each cell in the database. 
%
% function dataStruct=scc.diagnostic.buildImages(X)
%
%
% Purpose
% Helper function for buildImageDatabase.
% Call from project root dir to make the webpage. 
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

targetDir='./DiagnosticPlots/imageDatabase';
imageDir=fullfile(targetDir,'images');
if ~exist(imageDir,'dir')
    fprintf('Making directory %s\n', imageDir);
    mkdir(imageDir)
end


[visualAreaNames,colorMap]=brainAreaNames.visualAreas;

%First build top-view with the location of each neuron 
clf
set(gcf,'PaperPosition',[0,0,30,30])
aratools.projectAtlas.simplePlotter(projections(1))
hold on
for ii=1:length(X.data)
    sp=X.data(ii).pointsInARA.rawSparseData.sparsePointMatrix(1,:);
    plot(sp(2),sp(3),'*b')
    text(sp(2)+1,sp(3),num2str(ii),'color','r')
end
xlim([260,400])
ylim([90,238]) 

hold off
drawnow
somaFname=fullfile(imageDir,'soma.png');
print(gcf,'-dpng',somaFname)
fprintf('Saved soma image to %s\n',somaFname)


clf
set(gcf,'Color','w','Renderer','zbuffer')

% Loop through and make each cell plot
H=overlayCellsOnThreeProjections(X.data(1).pointsFname,projections,colorMap);

%Overlay sub-cortical areas present on this neuron
S=overlaySelectAreas(H,X.data(1).pointsInARA.leaves.ARAindex);


thisFname = fullfile(imageDir,X.data(1).details.cellID);
thisFname=[thisFname,'.png'];
print(gcf,'-dpng',thisFname)
fprintf('Saved image to %s\n',thisFname)

for ii=1:length(X.data)

    fprintf('%d/%d. %s\n', ii, length(X.data), X.data(ii).pointsFname)
    thisFname = fullfile(imageDir,X.data(ii).details.cellID);
    delete(S)
    H=overlayCellsOnThreeProjections(X.data(ii).pointsFname,H,colorMap);
    S=overlaySelectAreas(H,X.data(ii).pointsInARA.leaves.ARAindex);

    thisFname=[thisFname,'.png'];
    set(gcf,'PaperPosition',[0,0,30,29])
    print(gcf,'-dpng',thisFname)
    fprintf('Saved image to %s\n',thisFname)

end


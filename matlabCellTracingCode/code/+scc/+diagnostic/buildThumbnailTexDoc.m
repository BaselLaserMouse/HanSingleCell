function dataStruct=buildThumbnailTexDoc(X)
% Make a latex document that will contain all the V1 traced cells
%
% function dataStruct=scc.diagnostic.buildThumbnailTexDoc(X)
%
%
% Purpose
% Call from project root dir to make the latex document
% 



if nargin<1
    help(['scc.diagnostic.',mfilename])
    return
end






%Write the tex doc
absPath=fileparts(which('scc.diagnostic.buildThumbnailTexDoc'));
texHeaderPath = fullfile(absPath,'private/headerThumbs.tex');

targetDir='./DiagnosticPlots/thumbnailsDatabase';

if ~exist(targetDir,'dir')
    error('Expected to find %s directory. It is not there\n',targetDir)
end







texTargetPath=fullfile(targetDir,'v1_thumbnails_all.tex');
copyfile(texHeaderPath,texTargetPath)
doCells(X,texTargetPath,[],20,1)

OUT=scc.diagnostic.summariseTerminationTypes(X);

texTargetPath=fullfile(targetDir,'v1_thumbnails_clean.tex');
copyfile(texHeaderPath,texTargetPath)
doCells(X,texTargetPath, OUT(4).indexValuesOfcleanCells, 16, 1)


% Back-labelled cells only
BL=[];
for ii=1:length(X.data)
    if X.data(ii).details.isBackLabeled 
        BL(end+1)=ii;
    end
end
fprintf('\nFound %d Back-labelled cells\n', length(BL))
texTargetPath=fullfile(targetDir,'v1_backLabelled.tex');
if exist(texTargetPath,'file')
    delete(texTargetPath)
end
copyfile(texHeaderPath,texTargetPath)
doCells(X,texTargetPath,BL,length(BL),1)



% Non-V1 cells only
V1_IDs = [name2structureID('Primary visual area, layer 2/3'), ...
        name2structureID('Primary visual area, layer 1'), ... %Just in case the registration was a bit off. This is v. rare.
        name2structureID('Primary visual area, layer 4'), ...
        name2structureID('Primary visual area')];

noV1=OUT(5).indexValuesOfcleanCells;

fprintf('\nFound %d non-V1 cells\n', length(noV1))

texTargetPath=fullfile(targetDir,'non_v1_cells.tex');
if exist(texTargetPath,'file')
    delete(texTargetPath)
end
copyfile(texHeaderPath,texTargetPath)
doCells(X,texTargetPath,noV1,length(noV1),0)




function doCells(X,texTargetPath,IDs,thumbsPerFig,keepV1only)

    targetDir='./DiagnosticPlots/thumbnailsDatabase';
    imagePath=fullfile(targetDir,'images');
    d=dir(fullfile(imagePath,'*.eps'));
    fid = fopen(texTargetPath,'a+');



    %V1 cells only
    V1_IDs = [name2structureID('Primary visual area, layer 2/3'), ...
            name2structureID('Primary visual area, layer 1'), ... %Just in case the registration was a bit off. This is v. rare.
            name2structureID('Primary visual area, layer 4')];

    n=0; %Current number of panels in current figure
    cellNames = ''; %this contains the caption string
    for ii=1:length(X.data)
        %Skip excluded or non-V1 cells
        if  X.data(ii).details.excludeFromAnalysis
            continue
        end

        if ~isempty(IDs) && ~any(IDs==ii)
            continue
        end

        if keepV1only %If we are to exlcude cells outside of V1
            if X.data(ii).details.forceNotV1==1 || ~any(X.data(ii).pointsInARA.rawSparseData.ARAindex(1)==V1_IDs)
                continue
            end
        end


        thisEPS = fullfile(imagePath, [X.data(ii).details.cellID,'.eps']);


        if ~exist(thisEPS,'file')
            fprintf('Unable to find %s. SKIPPING\n', thisEPS)
            continue
        end

        fprintf('Doing cell %d/%d\n', ii, length(X.data))

        thisEPS = regexprep(thisEPS,'.*(images.*eps)','./$1');

        %Write figures
        if n==0
            fprintf(fid,'\n\\begin{figure}\n \\begin{adjustwidth}{-\\oddsidemargin-0.25in}{-\\rightmargin}\n \\begin{center}\n');
        end

        if thumbsPerFig>=20
            fprintf(fid,' \\includegraphics[width=1.25in]{%s}\n',thisEPS);
        else
            fprintf(fid,' \\includegraphics[width=1.5in]{%s}\n',thisEPS);
        end
        cellNames =[cellNames, X.data(ii).details.cellID, ', '];
        n=n+1;

        if n==thumbsPerFig || ii==length(X.data)
            cellNames(end-1:end)=[];
            cellNames = strrep(cellNames,'_','\_');
            fprintf(fid, ' \\end{center}\n \\end{adjustwidth}\n \\caption{%s}\n \\end{figure}\n\\clearpage\n',cellNames);
            n=0;
            cellNames='';
        end

    end

    writeFooterAndClose(fid)


% Convenience functions
function writeFooterAndClose(fid)
    fprintf(fid,'\n\n\\end{document}\n');
    fclose(fid);

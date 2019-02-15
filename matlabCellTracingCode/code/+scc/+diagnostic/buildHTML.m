function dataStruct=buildHTML(X)
% Make HTML pages only for the traced cells
%
% function dataStruct=scc.diagnostic.buildHTML(X)
%
%
% Purpose
% Call from project root dir to make the webpage
% helper function for buildImageDatabase
% 



if nargin<1
    help(['scc.diagnostic.',mfilename])
    return
end




targetDir='./DiagnosticPlots/imageDatabase';

if ~exist(targetDir,'dir')
    error('Expected to find %s directory. It is not there\n',targetDir)
end



%Write HTML
fid=openAndWriteHeader(fullfile(targetDir,'index.html'));

fprintf(fid,'<h2 id="top">Traced cell image database</h2>\n');
fprintf(fid,'<a href="%s"><img src="%s" width=800px /></a>\n\n\n', fullfile('./images','soma.png'), fullfile('./images','soma.png'));

%fprintf(fid,'<p>\n<a href="v1.html">V1 L2/3 cells only</a><br />\n\n');
%fprintf(fid,'<a href="non_v1.html">Non-V1 L2/3 cells only</a>\n</p>\n<hr />\n\n');



% Make all images
data = X.returnData;
cellList=scc.diagnostic.summariseTerminationTypes(X,true);

for ii=1:length(data)
    cellID = data(ii).details.cellID;
    if  data(ii).details.excludeFromAnalysis
        fprintf('Skipping cell %s -- marked for exlcusion\n', cellID)
        continue
    end
    thisPath=fullfile('./images',[cellID,'.png']);

    if any(cellList(1).indexValuesOfcleanCells == ii)
        cleanCells = 'No';
    else
        cleanCells = 'Yes';
    end

    if any(cellList(2).indexValuesOfAllCells == ii)
        inV1 = 'Yes';
    else
        inV1 = 'No';
    end

    if data(ii).details.isBackLabeled
        BL = 'Yes';
    else
        BL = 'No';
    end

    fprintf(fid, ...
        ['<h3>%d. %s</h3>\n', ...
         'Presence of premature terminations: %s \n<br />\n', ...
         'Soma in V1: %s \n<br />\n', ...
         'Back-labelled cell: %s \n<br />\n', ...
         '<a href="#top">Go to top</a></p>\n<a href="%s"><img src="%s" width=800px /></a>\n'],...
            ii, cellID, cleanCells, inV1, BL, thisPath, thisPath);
end

writeFooterAndClose(fid)


return
%Write only V1 cells in L2/3 (plus a little wiggle room in case a soma was mis-IDed as L1 or L4)
fid=openAndWriteHeader(fullfile(targetDir,'v1.html'));

fprintf(fid,'<h2 id="top">Traced cell image database: V1 only</h2>\n');



V1_IDs = [name2structureID('Primary visual area, layer 2/3'), ...
        name2structureID('Primary visual area, layer 1'), ... %Just in case the registration was a bit off. This is v. rare.
        name2structureID('Primary visual area, layer 4')];


%Go through each cell: and plot V1 only
for ii=1:length(C(2).indexValuesOfAllCells)
    if  X.data(ii).details.excludeFromAnalysis
        continue
    end
    if X.data(ii).details.forceNotV1==1 || ~any(X.data(ii).pointsInARA.rawSparseData.ARAindex(1)==V1_IDs)
        continue
    end

    thisPath=fullfile('./images',[X.data(ii).details.cellID,'.png']);
    fprintf(fid,'<h3>%d. %s</h3><a href="#top">Go to top</a></p>\n<a href="%s"><img src="%s" width=800px /></a>\n',...
            ii,X.data(ii).details.cellID, thisPath, thisPath);
end

writeFooterAndClose(fid)




%Write only non-V1 cells in L2/3 (plus a little wiggle room in case a soma was mis-IDed as L1 or L4)
fid=openAndWriteHeader(fullfile(targetDir,'non_v1.html'));

fprintf(fid,'<h2 id="top">Traced cell image database: non-V1 only</h2>\n');
fprintf('\n\nBuilding non-V1 cell page:\n')
%Go through each cell
n=0;
for ii=1:length(X.data)
    if  X.data(ii).details.excludeFromAnalysis
        continue
    end
    if X.data(ii).details.forceNotV1==0 && any(X.data(ii).pointsInARA.rawSparseData.ARAindex(1)==V1_IDs)
        continue
    end
    fprintf('%s\n',X.data(ii).details.cellID)
    thisPath=fullfile('./images',[X.data(ii).details.cellID,'.png']);
    fprintf(fid,'<h3>%d. %s</h3><a href="#top">Go to top</a></p>\n<a href="%s"><img src="%s" width=800px /></a>\n',...
            ii,X.data(ii).details.cellID, thisPath, thisPath);
    n=n+1;
end
fprintf('Found %d non-V1 cells\n',n)
writeFooterAndClose(fid)




% Convenience functions
function fid=openAndWriteHeader(fname)
    fprintf('Opening %s for writing\n', fname)
    fid = fopen(fname,'w+');
    if fid<0
        fprintf('Failed to open %s\n',fname)
        return
    end
    fprintf(fid,'<html>\n\n<head></head>\n\n\n<body>\n');

function writeFooterAndClose(fid)
    fprintf(fid,'\n<p><i>Made by scc.diagnostic.buildImageDatabase at %s</i></p>\n', datestr(now,'dd/mm/yyyy - HH:MM') );
    fprintf(fid,'</body>\n\n');
    fclose(fid);

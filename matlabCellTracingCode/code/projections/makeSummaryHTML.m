function makeSummaryHTML(copyFilesToOverlayDir)
% make a summary HTML file in the root dir of all experiments
% 
% function makeSummaryHTML(copyFilesToOverlayDir)
%
%
% Purpose
% Use the analysis_log.csv file to locate all processed animals
%
% expects an analysis_log.csv file:
% id,anim name,use/not use,notes
%
% 
% Inputs
% copyFilesToOverlayDir - [optional 0 by default] if 1, copy all png files to an overlay directory

if nargin<1
	copyFilesToOverlayDir=0;
end
if ~exist('analysis_log.csv')
	fprintf('Can not find analysis_log.csv\nCall this function from the root dir of all animals\n')
	return
end




%Get the directory where the overlays are located
S=settings_handler('settingsFiles_ARAtools.yml');
overlaysDir = fullfile(S.downSampledDir, S.sample2araDir, 'overlays');
dirs = aratools.utils.returnProcessedExperiments(overlaysDir);


%write the html file (and copy files if appropriate)
fid = fopen('cell_projections.html','w');
fprintf(fid,'<html>\n<head></head>\n<body>\n<table width=1000px border=1px>');

if copyFilesToOverlayDir
	copyDir = 'allOverlays';
	if ~exist(copyDir,'dir')
		mkdir(copyDir)
	end
end


n=1;
for ii=1:length(dirs)
	%get the names of the png files
	png=dir(fullfile(dirs{ii},'*.png'));
	if isempty(png)
		fprintf('no png files in %s\n',dirs{ii});
		continue
	end

	for jj=1:length(png)
		if mod(jj,3)==1
			fprintf(fid,'<tr><td><b>%d. %s</b></td>\n',n,strtok(dirs{ii},filesep));
			n=n+1;
		end

		pth2file=fullfile(dirs{ii},png(jj).name);
		fprintf(fid,'<td><a href="%s"><img src="%s" width=300px/></a></td>',pth2file,pth2file);

		if mod(jj,3)==0
			fprintf(fid,'</tr>\n');
		end

		if copyFilesToOverlayDir
			copyfile(pth2file,copyDir)
		end

	end

end

fprintf(fid,'</body>\n</table></html>');
function out=outlinesFromTemplate(template,doPlot)
% return coronal, sagittal, and transverse sections from the ARA template
%
% function out=outlinesFromTemplate(template,doPlot)	
%
% template is the matrix containing the template or path to template (mhd)
% doPlot shows the results and is 0 by default
%
%



if isstr(template)
	if exist(template,'file')
		template = mhd_read(template);
	end
end

if nargin<2
	doPlot=0;
end


transverse=(squeeze(max(template,[],1)));
sagittal=(squeeze(max(template,[],2)));
coronal=(squeeze(max(template,[],3)));

sections = {transverse, sagittal, coronal};


for ii=1:3
	bw{ii}=convertBW(sections{ii});
	bT{ii}=bwboundaries(bw{ii},'noholes');
end


if doPlot
	clf
	for ii=1:3
		subplot(2,2,ii)
		imagesc(bw{ii})
		hold on 
		plot(bT{ii}{1}(:,2),bT{ii}{1}(:,1),'-g','linewidth',5)
		hold off
	end
	colormap gray
end


out.sections = sections;
out.bw=bw;
out.traces=bT;


function bw=convertBW(im)
	bw = medfilt2(im,[2,2]);
	thresh=20;
	bw(bw<thresh)=0;
	bw(bw>=thresh)=1;


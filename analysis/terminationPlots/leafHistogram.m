function out=leafHistogram(data,doPlot)
% Distance of each leaf from the soma with abrupt terminations highlighted
%
%
%    function leafHistogram(data)
%
% Inputs
% data - one neuron from cleanCells.returnData object 
% doPlot - true by default
%
% e.g.
% >> load ~/tvtoucan/Mrsic-Flogel/hanyu/Analyses/cleanCells.mat
% >> D = cleanCells.returnData;
%
% ONE:
% >> OUT = leafHistogram(D(12))
%
% TWO: 
% for ii=1:length(D), OUT(ii) = leafHistogram(D(ii),false); end 
%
%
% Rob Campbell - Basel 2017
%
% See also:
% plotLeafData, plotLeafDataScatter



if nargin<2
    doPlot=true;
end


L=data.origTrace.findleaves;
voxelSize = str2num(data.voxelSize);

% Distance to root for each node
for ii=1:length(L)
    pth = data.origTrace.pathtoroot(L(ii));
    out.distToRoot(ii) = pathDistance(data.origTrace,pth,voxelSize); %distance to root
    
    % Which are premature terminations? (skip if no "nodeType" exists)
    if ~isfield(data.origTrace.Node{L(ii)}.data,'nodeType')
        out.premature(ii)=0;
        out.nodeType{ii}='normal';
        continue
    end

    nType = data.origTrace.Node{L(ii)}.data.nodeType;    

    if strcmp(nType,'normal')
        out.premature(ii)=0;
        out.nodeType{ii}='normal';
    else
        out.premature(ii)=1;
        out.nodeType{ii}=nType;
    end



end

% Number of branch nodes 
for ii=1:length(L)
    pth = data.origTrace.pathtoroot(L(ii)); % All nodes going back to the root node (the soma)
    for jj=length(pth):-1:1
        if length(data.origTrace.getchildren(pth(jj))) < 2
            pth(jj)=[];
        end
    end
    out.branchesToRoot(ii) = length(pth);
end



if doPlot
    cla
    hist(out.distToRoot,25)
    ptch=findobj(gca,'type','patch');
    set(ptch,'FaceColor',[1,1,1]*0.5)

    hold on
    for ii=1:length(out.premature)
        if out.premature(ii) ~= 0
            d=out.distToRoot(ii);
            plot([d,d],[0,2],'r-','linewidth',4)
        end
    end
    hold off

    [a,fn]=fileparts(data.pointsFname);

    str = sprintf('%s - %d premature terminations', fn, sum(find(out.premature>0)));
    str = regexprep(str,'_', '\\_');
    title(str)
end

out.hasPremature = any(out.premature);
out.cellID = data.details.cellID;

function d=pathDistance(treeData,pth,voxelSize)
    %Sum the euclidean distance along the path of points defined by pth
    d=0;

    for ii=1:length(pth)-1
        a = treeData.Node{pth(ii)};
        b = treeData.Node{pth(1+ii)};
        % Values are in microns apart from plane number, which is an index, and our planes are 10 microns apart. 
        m = [ [a.xVoxel, a.yVoxel, a.zVoxel*10] ; ...
              [b.xVoxel, b.yVoxel, b.zVoxel*10] ];
        d=d+pdist(m);
    end
    d=round(d);

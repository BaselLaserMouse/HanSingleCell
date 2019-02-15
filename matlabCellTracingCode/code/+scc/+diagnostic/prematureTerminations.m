function varargout=prematureTerminations(data,removeV1)
% function prematureTerminations(data,removeV1)
% 
% e.g.
% load sonastv/Data/Mrsic-Flogel/hanyu/Analyses/leavesByArea.mat
% prematureTerminations(data)
%
%
% Can we identify premature terminations in some easy way?
%
% TODO: appears not to work on leavesByArea -- there is no "distToRoot"

if nargin<2
    removeV1=0;
end

%Get plot data
leafDistR=[data.distToRoot];
leafDistN=[data.distToNearestBranch];
isWhite=[data.isWhiteMatter];
ind=vertcat(data.ind);

%We now determine for each sample, how many times each area appears
for ii=1:length(data)
    thisInd = data(ii).ind;
    data(ii).counts=zeros(size(thisInd));
    u=unique(data(ii).ind)';
    for k=u
        f=thisInd==k;
        data(ii).counts(f) = sum(f);
    end
end
counts=vertcat(data.counts);

if removeV1
    labels=getAllenStructureList;
    f=strmatch('Primary visual',labels.name);
    VISpInd = labels.id(f);

    f=ismember(ind,VISpInd);
    fprintf('Removing %d V1 leaves\n',sum(f))
    
    leafDistR(f)=[];
    leafDistN(f)=[];
    isWhite(f)=[];
    counts(f)=[];
    excludeText=' (excluding V1)';
else
    excludeText='';
end




clf

subplot(2,2,1)
plot(leafDistR,leafDistN,'.k')
hold on 
plot(leafDistR(isWhite),leafDistN(isWhite),'or')
hold off
L=unityLine;
L.Color='g';
title(sprintf('%d leaves from %d samples%s',length(leafDistN),length(data),excludeText))
xlabel('Leaf distance from soma')
ylabel('Leaf distance to nearest branch')


subplot(2,2,2)
plot(counts,leafDistN,'.k')
hold on 
plot(counts(isWhite),leafDistN(isWhite),'or')
hold off
title(sprintf('%d leaves from %d samples%s',length(leafDistN),length(data),excludeText))
xlabel('Counts in area')
ylabel('Leaf distance to nearest branch')


subplot(2,2,3)
plot3(leafDistR,leafDistN,counts,'.k')
hold on 
plot3(leafDistR(isWhite),leafDistN(isWhite),counts(isWhite),'or')
hold off
xlabel('Leaf distance from soma')
ylabel('Leaf distance to nearest branch')
zlabel('Counts in area')
grid on


if nargout>0
    varargout{1}=data;
end
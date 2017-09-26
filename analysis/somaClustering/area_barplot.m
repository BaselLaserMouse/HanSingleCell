function varargout = area_barplot(data,errB)
    % Bar graph showing how many cells project to each area; subdivide by ML and RL
    % 
    % function varargout = area_barplot(data)
    %
    % Where "data" is:
    %
    % >> for ii=2:13, [~,data{ii}]=clusterPos(allCellMat,cleanCells,ii); end
    %
    % errB is the bootstrapped error bar map produced by this function.
    % feeding it in as an input argument stops it being regenerated each time.
    %
    % Rob Campbell - Basel 2017

data(cellfun(@isempty,data))=[];

%Minus 2 because some brain areas get no connections
pltData.all = zeros(1,length(data)-2);
pltData.R = zeros(1,length(data)-2);
pltData.C = zeros(1,length(data)-2);
pltData.M = zeros(1,length(data)-2);
pltData.L = zeros(1,length(data)-2);

% to help make the bar labels 
[n,c,abrv]=brainAreaNames.visualAreas;
xLabels = containers.Map(c.areaNames, abrv(1:end-1) );

n=1;
for ii=1:length(data)
    tD = data{ii};

    f = ~cellfun(@isempty,tD.Target);

    tmpInd = find(f);
    if isempty(tmpInd)
        fprintf('Skipping %d\n',ii)
        continue
    end
    ind = tmpInd(1);
    pltData.targetName{n} = tD.Target{ind};
    pltData.abrvName{n} = xLabels(tD.Target{tmpInd(1)});

    pltData.all(n) = sum(f);

    pltData.R(n) = sum(f & cellfun(@(x) strcmp('R',x),tD.RC)); %To this area and rostral
    pltData.C(n) = sum(f & cellfun(@(x) strcmp('C',x),tD.RC)); %To this area and caudal
    pltData.M(n) = sum(f & cellfun(@(x) strcmp('M',x),tD.ML)); 
    pltData.L(n) = sum(f & cellfun(@(x) strcmp('L',x),tD.ML));

    n=n+1;
end

% Let's bootstrap some error bars based on the sample sizes
if nargin<2
    n=unique(pltData.all);
    n(n<4)=[];
    errB = containers.Map('KeyType','double','ValueType','any');
    for ii=1:length(n)
        errB(n(ii))=btstrpErrBar(data{1},n(ii),10E3); % generate a bootstrap chance value for this sample size
    end
end


%So the labels for the x axis are:
fSize=15;
clf 

subplot(2,1,1)
title('R/C')
hR=bar(pltData.R);
hold on 
hC=bar(-pltData.C);
ylim([-20,20])

set(hR,'FaceColor',[1,0.55,0.55],'EdgeColor',[0.75,0,0])
set(hC,'FaceColor',[0.55,0.55,1],'EdgeColor',[0,0,0.75])

W = 0.35; % line width for error bars
erbMrkr = {'k:','LineWidth',2};
% Add the error bars
for x = 1:length(pltData.all)
    n = pltData.all(x); %total points here
    if n<4, continue, end

    mu=mean(errB(n).R);
    plot([x-W,x+W], [mu,mu], erbMrkr{:})
    mu=-mean(errB(n).C);
    S=std(errB(n).C)*1.96;
    plot([x-W,x+W], [mu,mu], erbMrkr{:})


    %Look for significant differences
    obs=abs(pltData.R(x)-pltData.C(x));
    bs = errB(n).R-errB(n).C;
    p=length(find(bs>obs))/length(bs);

    if p<0.025
        text(x-0.4, pltData.R(x)+2,sprintf('p=%0.4f',p),'FontSize',11)
    end

end
set(gca,'XTickLabel',pltData.abrvName,'YTick',-20:10:20,...
        'YTickLabel',[20,10,0,10,20],'FontSize',fSize)
box off
text(-0.65,10,'Rostral', 'Rotation',90,'FontSize',fSize)
text(-0.65,-10,'Caudal', 'Rotation',90,'FontSize',fSize)



subplot(2,1,2)
title('M/L')
hM=bar(pltData.M);
hold on 
hL=bar(-pltData.L);
ylim([-20,20])

set(hM,'FaceColor',[1,0.55,0.55],'EdgeColor',[0.75,0,0])
set(hL,'FaceColor',[0.55,0.55,1],'EdgeColor',[0,0,0.75])


% Add the error bars
for x = 1:length(pltData.all)
    n = pltData.all(x); %total points here
    if n<4, continue, end

    mu=mean(errB(n).M);
    plot([x-W,x+W], [mu,mu], erbMrkr{:})
    mu=-mean(errB(n).L);
    plot([x-W,x+W], [mu,mu], erbMrkr{:})

    %Look for significant differences
    obs=abs(pltData.M(x)-pltData.L(x));
    bs = errB(n).M-errB(n).L;
    p=length(find(bs>obs))/length(bs);

    if p<0.025
        text(x-0.4,pltData.M(x)+2,sprintf('p=%0.4f',p),'FontSize',11)
    end

end
set(gca,'XTickLabel',pltData.abrvName,'YTick',-20:10:20,...
    'YTickLabel',[20,10,0,10,20],'FontSize',fSize)
box off
text(-0.65,10,'Medial', 'Rotation',90,'FontSize',fSize)
text(-0.65,-10,'Lateral', 'Rotation',90,'FontSize',fSize)

if nargout>0
    varargout{1}=pltData;
end
if nargout>1
    varargout{2}=errB;
end


function errB = btstrpErrBar(d,n,reps) % generate a bootstrap chance value for this sample size
    % errB(1) is the number on the rostral side
    % errB(2) is the number on the ML

    errB.R = zeros(1,reps);
    errB.C = zeros(1,reps);
    errB.M = zeros(1,reps);
    errB.L = zeros(1,reps);

    for ii=1:reps
        r = randperm(length(d.CellID));
        r = r(1:n);
        errB.R(ii) = sum( cellfun(@(x) strcmp('R',x), d.RC(r)) ); %To this area and rostral
        errB.C(ii) = sum( cellfun(@(x) strcmp('C',x), d.RC(r)) ); %To this area and caudal
        errB.M(ii) = sum( cellfun(@(x) strcmp('M',x), d.ML(r)) ); 
        errB.L(ii) = sum( cellfun(@(x) strcmp('L',x), d.ML(r)) ); 
    end


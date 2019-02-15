function pProjectionsByArea(pltOut)
% Plot the probability of finding a projection in a given brain area
%
%
% Inputs
% pltOut - the output of pointsByAreaPlot. So this function will only process
%          areas that were plotted. If you exlcuded areas and they weren't on
%          the plot, then they won't be shown here. 
%
%
% Also see:
% brainAreaBarChart



% Plot the number of projections in each brain area. So if a cell terminates 10 
% times in a particular area it is counted 10 times. 
clf

ax1=axes('position',[0.1,0.54,0.8,0.3]);
data=sum(pltOut.dataMat,2);
[~,ind]=sort(data);

plot(data(ind),'ob-','MarkerFaceColor',[0.5,0.5,1])
set(gca,'XTick',1:length(data),'XTickLabel',[])
grid on
ylabel('Number of projections')

%Overlay the probability
overlayProbability(ax1,data)


%Add little histogram showing the number of projection targets per cell
axes('position',[0.125,0.675,0.3,0.13])

%binarize
b=pltOut.dataMat;
b(b>0)=1;

data=sum(b,1);
hist(data)
set(gca,'YAxisLocation','right')
xlabel('# proj. targets')

% Plot the number of projections in each brain area. So if a cell terminates 10 
% times in a particular area it is counted 10 times. 
ax2=axes('position',[0.1,0.21,0.8,0.3]);
data=sum(b,2); 

%sort the same way as above
plot(data(ind),'ob-','MarkerFaceColor',[0.5,0.5,1])
set(gca,'XTick',1:length(data),'XTickLabel',pltOut.areaNamesInSamples(ind),'XTickLabelRotation',45)
grid on
ylabel('Number of neurons')



%Overlay the probability
overlayProbability(ax2,data)




%--------------------------------------
function overlayProbability(ax,data)
    axes('Position',get(ax,'position'))
    set(gca,'YAxisLocation','right','Color','None','XTick',[])

    yIN = get(ax,'YLim');
    Y=round([0,yIN(2)/sum(data(:))],2);

    ylim(Y);
    
    nYTicks=length(get(ax,'YTick'));
    set(gca,'YTick',linspace(Y(1),Y(2),nYTicks))
    YL=get(gca,'YTickLabel');

    set(gca,'YTickLabel', cellfun(@(x) round(str2double(x),2), YL) )
    ylabel('probability')

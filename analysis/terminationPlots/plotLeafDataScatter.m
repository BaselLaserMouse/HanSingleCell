function out=plotLeafDataScatter(data)
% Plot output of plotLeafDataScatter
%
%    function plotLeafData(data)
%
% data - the output of leafHistogram run over all cells as follows:
%        for ii=1:length(D), OUT(ii) = leafHistogram(D(ii),false); end 
%
%
% Rob Campbell - Basel 2017
%
% See also:
% leafHistogram, plotLeafData

mrkr = {'ok', 'markerfacecolor',[1,1,1]*0.5, 'markersize',7};

clf


% Plot distance from soma of each premature termination as a function of the furthest non-premature distance from the soma. 

subplot(1,2,1)
hold on

pltData=[];

n=0;
for ii=1:length(data)
    indPremature = find(data(ii).premature);

    if isempty(indPremature)
        continue 
    end

    x = data(ii).distToRoot(find(~data(ii).premature));

    for jj=1:length(indPremature)
        switch data(ii).nodeType{indPremature(jj)}
        case 'fading'
            fprintf('Fading termination in cell %s\n', data(ii).cellID)
            mrkr={'or','MarkerFaceColor',[1,0.7,0.7]};
        case 'callosal_fading'
            fprintf('Fading callosal termination in cell %s\n', data(ii).cellID)
            mrkr={'or','MarkerFaceColor',[1,0.2,0.2]};
        otherwise
            fprintf('%s termination in cell %s\n', data(ii).nodeType{indPremature},data(ii).cellID)
            mrkr={'ok'};
        end

        %% COMMENT OUT TO HIGHLIGHT FADING TERMINATIONS. 
        mrkr={'ok','markerfacecolor', [1,1,1]*0.5}; 

        y = data(ii).distToRoot(indPremature(jj));
        plot(max(x),y,mrkr{:})

        n=n+1;
        pltData(n,1) = y - max(x);
    end
end

title( sprintf('%d abrupt terminations',size(pltData,1)) )
X = xlim;

Y = ylim;
M=max([X(2),Y(2)]);
ylim([0,M])
xlim([0,M])
plot([0,M],[0,M],'--')
grid on
axis square
box on
xlabel('Max distance to root non-premature (\mum)')
ylabel('Max distance to root premature only (\mum))')


subplot(1,2,2)
hold on
n=0;
for ii=1:length(data)
    y = data(ii).branchesToRoot(find(data(ii).premature));
    if isempty(y)
        continue 
    end
    x = data(ii).branchesToRoot(find(~data(ii).premature));

    % Get number of branches to root for the furthest terminals
    [~,indX] = max( data(ii).distToRoot(find(~data(ii).premature)) );
    [~,indY] = max( data(ii).distToRoot(find(data(ii).premature)) );

    plot(x(indX), y(indY), mrkr{:})
    n=n+1;
    pltData(n,2) = median(y) - median(x);
end
title( sprintf('%d abrupt terminations',size(pltData,1)) )

X = xlim;
Y = ylim;
M=max([X(2),Y(2)]);
ylim([0,M])
xlim([0,M])
plot([0,M],[0,M],'--')
grid on
axis square
box on
xlabel('Number of branches to root (furthest non-premature)')
ylabel('Number of branches to root (furthest premature)')

jitter(gca,2,true)
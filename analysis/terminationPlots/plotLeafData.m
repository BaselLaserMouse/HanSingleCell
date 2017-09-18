function out=plotLeafData(data)
% Plot output of plotLeafData
%
% function plotLeafData(data)
% 
% data - the output of leafHistogram run over all cells as follows:
%        for ii=1:length(D), OUT(ii) = leafHistogram(D(ii),false); end 
%
%
% Rob Campbell - Basel 2017
%
% See also:
% leafHistogram, plotLeafDataScatter

histBins = 100;

d = [data.distToRoot]/1E3;

%Pull out the distances of the premature terminals
p = d(find([data.premature]));

%normalise
for ii=1:length(data)
    data(ii).distToRootNorm = data(ii).distToRoot/max(data(ii).distToRoot);
end
dn = [data.distToRootNorm];

pn = dn(find(p)); %premature


clf

subplot(2,4,1)
hist(d,histBins)
ptch=findobj(gca,'type','patch');
set(ptch,'FaceColor',[1,1,1]*0.5)



Y=ylim;
hold on
for ii=1:length(p)
    t = p(ii);
    plot([t,t],[0,Y(2)*0.05], '-r', 'LineWidth', 2)
end

xlabel('Distance from soma (mm)')
ylabel('n')


%normalised
subplot(2,4,2)


hist(dn,histBins)
ptch=findobj(gca,'type','patch');
set(ptch,'FaceColor',[1,1,1]*0.5)


Y=ylim;
hold on
for ii=1:length(pn)
    t = pn(ii);
    plot([t,t],[0,Y(2)*0.05], '-r', 'LineWidth', 2)
end

xlabel('Distance from soma (normalised')
ylabel('n')





% Cumulative plots of the above two
subplot(2,4,5)

% What proportion of terminals falls into each size?
sizes = 0:0.1:12;
n = zeros(1,length(sizes));

for ii = 1:length(sizes)
    n(ii) = length(find(d<=sizes(ii)) );
end
n = n / max(n);
plot(sizes,n ,'-k','linewidth',2)
set(gca,'XTick', 0:14)
ylim([0,max(n)*1.05])

%Y=ylim;
hold on

%same thing for the premature terminations
np = zeros(1,length(sizes));
for ii=1:length(sizes)
    np(ii) = length(find(p<=sizes(ii)) );
end
np = np / max(np);

plot(sizes,np,'-r','linewidth',2)

grid on
hold off

xlabel('Distance from soma (mm)')
ylabel('Cumulative proportion')




% Normalised
subplot(2,4,6)

% What proportion of terminals falls into each size?
sizes = 0:0.05:1;
n = zeros(1,length(sizes));

for ii = 1:length(sizes)
    n(ii) = length(find(dn<=sizes(ii)) );
end
n = n / max(n);
plot(sizes,n ,'-k','linewidth',2)
set(gca,'XTick', 0:0.1:1)
ylim([0,max(n)*1.05])

%Y=ylim;
hold on

%same thing for the premature terminations
npn = zeros(1,length(sizes));
for ii=1:length(sizes)
    npn(ii) = length(find(pn<=sizes(ii)) );
end
npn = npn / max(npn);

plot(sizes,npn,'-r','linewidth',2)

grid on
hold off

xlabel('Distance from soma (normalised)')
ylabel('Cumulative proportion')




%----------------------------------
% Plot the number of branches between each terminal node and the soma
b = [data.branchesToRoot];

%Pull out the values of the premature terminals
pb = b(find([data.premature]));

%normalise
for ii=1:length(data)
    data(ii).branchesToRootNorm = data(ii).branchesToRoot/max(data(ii).branchesToRoot);
end
bn = [data.branchesToRootNorm];

pbn = bn(find(p)); %premature normalised



subplot(2,4,3)
hist(b,10)
ptch=findobj(gca, 'type', 'patch');
set(ptch, 'FaceColor', [1,1,1]*0.5)



Y=ylim;
hold on
for ii=1:length(pb)
    t = pb(ii);
    plot([t,t],[0,Y(2)*0.05], '-r', 'LineWidth', 2)
end

xlabel('# Branches to soma')


% Normalised
subplot(2,4,4)
hist(bn,10)
ptch=findobj(gca, 'type', 'patch');
set(ptch, 'FaceColor', [1,1,1]*0.5)



Y=ylim;
hold on
for ii=1:length(pb)
    t = pbn(ii);
    plot([t,t],[0,Y(2)*0.05], '-r', 'LineWidth', 2)
end

xlabel('# Branches to soma')


% Cumulative plots of the above two
subplot(2,4,7)

% What proportion of terminals falls into each size?
sizes = 0:0.5:28;
n = zeros(1,length(sizes));

for ii = 1:length(sizes)
    n(ii) = length(find(b<=sizes(ii)) );
end
n = n / max(n);
plot(sizes,n ,'-k','linewidth',2)
set(gca,'XTick', 0:2:28)
ylim([0,max(n)*1.05])

%Y=ylim;
hold on

%same thing for the premature terminations
np = zeros(1,length(sizes));
for ii=1:length(sizes)
    np(ii) = length( find(pb<=sizes(ii)) );
end
np = np / max(np);

plot(sizes,np,'-r','linewidth',2)

grid on
hold off

xlabel('# Branches to soma')
ylabel('Cumulative proportion')



subplot(2,4,8)

% What proportion of terminals falls into each size?
sizes = 0:0.05:1;
n = zeros(1,length(sizes));

for ii = 1:length(sizes)
    n(ii) = length(find(bn<=sizes(ii)) );
end
n = n / max(n);
plot(sizes,n ,'-k','linewidth',2)
set(gca,'XTick', 0:0.1:1)
ylim([0,max(n)*1.05])

%Y=ylim;
hold on

%same thing for the premature terminations
np = zeros(1,length(sizes));
for ii=1:length(sizes)
    np(ii) = length( find(pbn<=sizes(ii)) );
end
np = np / max(np);

plot(sizes,np,'-r','linewidth',2)

grid on
hold off

xlabel('# Branches to soma (normalised)')
ylabel('Cumulative proportion')




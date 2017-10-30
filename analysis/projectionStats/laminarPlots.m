function varargout = laminarPlots(laminarData,cleanCells)
    % Plot distribution of terminals by layer within a subset of defined areas
    %
    % Distribution plots
    % 
    % out = laminarPlots(laminarData,cleanCells)
    %
    %
    % Inputs
    % cleanCells - xylem structure
    % laminarData - produced by getLaminarData
    %

    cleanOnly=false;
    if cleanOnly
        % Filter neurons by premature/not 
        s=scc.diagnostic.summariseTerminationTypes(cleanCells,true);
        s=s(4); %No back-labelled cells V1 only
        cleanCells = s.IDsofCleanCells;

        %Go through everything and remove non-clean cells
        for ii=1:length(laminarData)
            tData=laminarData(ii);
            ind=[];
            for jj=1:length(tData.cellID)
                if ~isempty( strmatch(tData.cellID{jj},cleanCells) )
                    ind(end+1)=jj;
                end
            end % jj
            tData.cellID = tData.cellID(ind);
            tData.counts = tData.counts(:,:,ind);
            laminarData(ii) = tData;
        end
    end


    disp('CHECK UNITS ON Y AXIS')



    doBoxPlots
    %doLinePlots

    labelEdgeSubPlots('Layer',sprintf('Proportion of axon length in each layer\nnormalized by layer cross-section area') )

    c=get(gcf,'children');
    yl=[];
    for ii=length(c):-1:1
        if ~strcmp(c(ii).Type,'axes')
            c(ii)=[];
            continue
        end
        tmp = c(ii).YLim;
        yl = [yl;tmp];
    end

    set(c,'YLim',[0, max(yl(:,2))])
    set(gcf,'PaperPosition',[0,0,40,30])





    function doLinePlots
        % Linked points as lines so we know which cell is which in all layers
        clf
        n=1;
        for ii=1:length(laminarData)


            if isempty(laminarData(ii).counts)
                fprintf('Nothing for area %s\n', laminarData(ii).areaName)
                continue
            end
            if size(laminarData(ii).counts,3)<3
                fprintf('Skipping  area %s with only %d cells\n', laminarData(ii).areaName, ...
                    size(laminarData(ii).counts,3)<3)
                continue
            end

            subplot(3,4,n)

            tmp = squeeze(laminarData(ii).counts(:,3,:));
            tmp = single(tmp);
            tmp = tmp ./ sum(tmp,1);
            [~,ind]=sort(max(tmp,[],1),'descend');

            plot(tmp)
            set(gca,'defaultAxesColorOrder',parula(length(ind)))
%            title( sprintf('%s (%d cells)',...
%                laminarData(ii).areaName, size(laminarData(ii).counts,3)) )
            title( sprintf(laminarData(ii).shortName))

            % Make nice x axis labels
            lab=cellfun(@(x) regexprep(x,'.* ',''), laminarData(ii).areaTable.name ,'UniformOutput', false);
            set(gca,'XTickLabel',lab)
            n=n+1;
            grid on
            box off
        end

    end %doLinePlots




    function doBoxPlots(rowToPlot)
        % Make notBoxPlots. One per layer
        clf
        n=1;

        for ii=1:length(laminarData)


            if isempty(laminarData(ii).counts)
                fprintf('Nothing for area %s\n', laminarData(ii).areaName)
                continue
            end

            if size(laminarData(ii).counts,3)<10
                if size(laminarData(ii).counts,3)==1
                    plural='';
                else
                    plural='s';
                end

                fprintf('Skipping "%s" with only %d cell%s.\n', laminarData(ii).areaName, ...
                    size(laminarData(ii).counts,3)<3, plural)
                continue
            end

            %Remove any cells that need excluding, should this have been requested
            if isfield(laminarData(ii),'brainsToExclude') && ~isempty(laminarData(ii).brainsToExclude)
                excl = laminarData(ii).brainsToExclude;
                indToRemove=[];
                for ff=1:length(excl)
                    indToRemove=[indToRemove;strmatch(excl{ff},laminarData(ii).cellID)];
                end
                laminarData(ii).counts(:,:,indToRemove)=[];
                laminarData(ii).cellID(indToRemove)=[];
                laminarData(ii).traceFiles(indToRemove)=[];
            end

            subplot(2,4,n)
            counts = squeeze(laminarData(ii).counts(:,3,:))';
            counts = single(counts);

            % Because we will group layer 6a and 6ba
            if size(counts,2)~=6
                fprintf('Skipping "%s" where only %d layers have been extracted.\n', laminarData(ii).areaName, ...
                    size(counts,2))
                continue
            end

            counts(:,5) = counts(:,5) + counts(:,6); counts(:,6)=[];

            % turn into proportion


            vol = squeeze(laminarData(ii).counts(:,2,:))';
            vol = single(vol);
            vol(:,5) = vol(:,5) + vol(:,6); vol(:,6)=[];
            vol = vol/length(laminarData(ii).planesInARA); %2D layer area

            if 1
                plotData = counts./vol; %Calculate the density 
                plotData = plotData ./ sum(plotData,2);
            else
                counts=counts ./ sum(counts,2);
            end

            H=notBoxPlot(plotData);
            set([H.data],'MarkerSize',3,'MarkerFaceColor',[1,1,1]*0.25,'LineWidth',0.25)
            delete([H.sdPtch])
            set([H.semPtch], 'FaceColor', [0.75,0.67,0.67], 'EdgeColor', [0.66,0.5,0.5])
            set([H.mu],'LineWidth',3)
            %Add medians

            M=median(plotData);
            hold on

            for kk=1:length(M)
                y=[M(kk),M(kk)];
                x=H(kk).mu.XData;
                plot(H(kk).mu.XData,y,'-b','LineWidth',H(kk).mu.LineWidth)
            end

%            title( sprintf('%s (%d cells)',...
%                laminarData(ii).areaName, size(laminarData(ii).counts,3)) )
            title( sprintf(laminarData(ii).shortName))

            % Make nice x axis labels
            lab=cellfun(@(x) regexprep(x,'.* ',''), laminarData(ii).areaTable.name ,'UniformOutput', false);
            lab{5} = '6';
            set(gca,'XTickLabel',lab(1:5))
            n=n+1;
            %ylim([-100,4E3])
        end

    end %doBoxPlots

end %close laminarPlots
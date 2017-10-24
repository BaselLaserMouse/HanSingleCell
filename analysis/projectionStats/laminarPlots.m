function varargout = laminarPlots(laminarData,cleanCells,cleanOnly)
    % Plot distribution of terminals by layer within a subset of defined areas
    %
    % out = laminarPlots(laminarData,cleanCells,cleanOnly)
    %
    %
    % Inputs
    % laminarData - see below
    % cleanCells - xylem structure
    % cleanOnly - false by default

    if nargin<3
        cleanOnly=false;
    end

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


    
    labelEdgeSubPlots('Layer','axon length') 

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

            subplot(3,3,n)

            tmp = squeeze(laminarData(ii).counts(:,2,:));
            plot(tmp)

            title( sprintf('%s (%d cells)',...
                laminarData(ii).areaName, size(laminarData(ii).counts,3)) )

            % Make nice x axis labels
            lab=cellfun(@(x) regexprep(x,'.* ',''), laminarData(ii).areaTable.name ,'UniformOutput', false);
            set(gca,'XTickLabel',lab)
            n=n+1;
            grid on
            box off
        end

    end %doLinePlots




    function doBoxPlots
        % Make notBoxPlots. One per layer
        clf
        n=1;
        for ii=1:length(laminarData)


            if isempty(laminarData(ii).counts)
                fprintf('Nothing for area %s\n', laminarData(ii).areaName)
                continue
            end
            if size(laminarData(ii).counts,3)<3
                fprintf('Skipping  area %s with only %disp cells\n', laminarData(ii).areaName, ...
                    size(laminarData(ii).counts,3)<3)
                continue
            end

            subplot(3,3,n)
            counts = squeeze(laminarData(ii).counts(:,3,:))';
            counts = single(counts);

            vol = squeeze(laminarData(ii).counts(:,2,:))';
            vol = single(vol);
            vol = vol/length(laminarData(ii).planesInARA);

            counts = counts*5; % TODO - should be correct but check
            if ~isvector(counts)
                notBoxPlot(counts./vol)
            else
                plot(counts,'ok-')
            end
            title( sprintf('%s (%d cells)',...
                laminarData(ii).areaName, size(laminarData(ii).counts,3)) )

            % Make nice x axis labels
            lab=cellfun(@(x) regexprep(x,'.* ',''), laminarData(ii).areaTable.name ,'UniformOutput', false);
            set(gca,'XTickLabel',lab)
            n=n+1;
            %ylim([-100,4E3])
        end

    end %doBoxPlots

end %close laminarPlots
function varargout = laminarPlots(laminarData,cleanCells,cleanOnly)
    % Plot distribution of terminals by layer within a subset of defined areas
    %
    % out = laminarPlots(laminarData,cleanOnly)
    %
    %
    % Inputs
    % laminarData - see below
    % cleanCells - xylem structure
    % cleanOnly - false by default
    %
    % Example
    % >> for ii=1:length(visNames); OUT(ii)=getLaminarData(cleanCells,visNames{ii}, LOADED_ARA); end
    % >> laminarPlots(OUT)

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



    %doBoxPlots
    doLinePlots


    function doLinePlots
        clf
        n=1;
        for ii=1:length(laminarData)


            if isempty(laminarData(ii).counts)
                fprintf('Nothing for area %s\n', laminarData(ii).areaName)
                continue
            end
            if size(laminarData(ii).counts,3)<3
                fprintf('Skipping  area %s with only % cells\n', laminarData(ii).areaName, ...
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

        labelEdgeSubPlots('Layer','axon length') 
    end %doLinePlots


    function doBoxPlots
        clf
        n=1;
        for ii=1:length(laminarData)


            if isempty(laminarData(ii).counts)
                fprintf('Nothing for area %s\n', laminarData(ii).areaName)
                continue
            end
            if size(laminarData(ii).counts,3)<3
                fprintf('Skipping  area %s with only % cells\n', laminarData(ii).areaName, ...
                    size(laminarData(ii).counts,3)<3)
                continue
            end

            subplot(3,3,n)
            tmp = squeeze(laminarData(ii).counts(:,2,:))';
            tmp = single(tmp);
            if ~isvector(tmp)
                notBoxPlot(tmp)
            else
                plot(tmp,'ok-')
            end
            title( sprintf('%s (%d cells)',...
                laminarData(ii).areaName, size(laminarData(ii).counts,3)) )

            % Make nice x axis labels
            lab=cellfun(@(x) regexprep(x,'.* ',''), laminarData(ii).areaTable.name ,'UniformOutput', false);
            set(gca,'XTickLabel',lab)
            n=n+1;
        end

        labelEdgeSubPlots('Layer','axon length') 
    end %doBoxPlots

end %close laminarPlots
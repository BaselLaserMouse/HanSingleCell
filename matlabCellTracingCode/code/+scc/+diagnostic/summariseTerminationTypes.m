function out=summariseTerminationTypes(X,quiet)
    % Summarise the number of neurons exhibiting different sorts of
    % seemingly abrupt terminations. i.e. those that don't branch
    % where we might a priori expect them to do so. 
    %
    % function out = scc.diagnostic.summariseTerminationTypes(X,quiet)
    %
    % Inputs
    % X - the xylem data object
    % quiet - false by default
    %
    %
    % Example:
    % >> OUT = scc.diagnostic.summariseTerminationTypes(X);
    % The data index values for V1 neurons:
    % >> OUT(2).indexValuesOfcleanCells
    % So you can do:
    % >> d=X.returnData;
    % >> cleanCells = d(OUT(2).indexValuesOfcleanCells);
    %
    %
    % The cells used in the manuscript are: 
    % unique([out(1).indexValuesOfcleanCells,out(3).indexValuesOfAllCells]) 
    % i.e. all clean cells and all the striatum cells

    if nargin<2
        quiet=false;
    end

    data=X.returnData;

    out(1)=searchCells(data,0);   % All cells
    out(2)=searchCells(data,1);    % V1 only
    out(3)=searchCells(data,0,1); % Only backlabelled
    out(4)=searchCells(data,1,-1); % Not backlabelled V1 only
    out(5)=searchCells(data,-1,-1); % Not backlabelled non-V1 only

    if ~quiet
        for ii=1:length(out)
           prettyDisplay(out(ii))
        end
    end




    function out=searchCells(data,V1only,backlabelled)

        if nargin<2 || isempty(V1only)
            %if zero we let through both V1 and non-V1
            %if 1 only V1
            %if -1 only non-V1
            V1only=0;
        end

        if nargin<3 || isempty(backlabelled)
            %if zero we let through all cells
            %if 1 only backlabelled
            %if -1 not backlabelled
            backlabelled=0;
        end

        if backlabelled>0
            out.backlabelled='only backlabelled cells';
        elseif backlabelled<0
            out.backlabelled='No backlabelled cells';
        elseif backlabelled==0
            out.backlabelled='Both backlabelled and non-backlabelled cells';
        end

        numAbruptTerminations=zeros(1,length(data));
        numCallosal=zeros(1,length(data));
        numBright=zeros(1,length(data));
        numBrightCallosal=zeros(1,length(data));
        numFading=zeros(1,length(data));
        numWhiteMatter=zeros(1,length(data));
        numWhiteMatterNonCallosal=zeros(1,length(data));
        numFadingGrayMater=zeros(1,length(data));
        cleanCell=zeros(1,length(data)); %All cells with a 1 are clean cells with no abrupt terminations
        allCells=zeros(1,length(data)); %All cells with a 1 are those that meet the selection criteria but may not be clean
        wInds=aratools.utils.whiteMatterInds;


        %In case we want to keep only V1 cells
        V1_IDs = [name2structureID('Primary visual area, layer 2/3'), ...
        name2structureID('Primary visual area, layer 1'), ... %Just in case the registration was a bit off. This is v. rare.
        name2structureID('Primary visual area, layer 4'), ...
        name2structureID('Primary visual area')];


        n=0;
        for ii=1:length(data)


            if V1only==1 && (data(ii).details.forceNotV1==1 || ~any(data(ii).pointsInARA.rawSparseData.ARAindex(1)==V1_IDs))
                continue
            elseif V1only==-1 && any(data(ii).pointsInARA.rawSparseData.ARAindex(1)==V1_IDs) && data(ii).details.forceNotV1==0
                continue
            end


            if backlabelled<0 && data(ii).details.isBackLabeled==true
                continue
            end

            if backlabelled>0 && data(ii).details.isBackLabeled==false
                continue
            end

            oTrace=data(ii).origTrace;
            n=n+1;

            %What sort of endings are there?
            for jj=1:length(oTrace.Node)

                if ~isfield(oTrace.Node{jj}.data,'nodeType')
                    %Skip if no node type exists for this node
                    continue
                end

                nType=oTrace.Node{jj}.data.nodeType;

                if strcmp(nType,'normal')
                    continue
                end

                %If we are here, this is a node which terminated abruptly
                numAbruptTerminations(ii) = numAbruptTerminations(ii)+1; %Count the number of abrupt terminations in this cell

                %What kind of termination was it? (bright or fading)
                if ~isempty(regexp(nType,'bright','once'))
                    numBright(ii)=numBright(ii)+1;
                end
                if ~isempty(regexp(nType,'fading','once'))
                    numFading(ii)=numFading(ii)+1;
                end

                %Count the number of white matter endings in the corpus callosum and outside of it
                if any(wInds==data(ii).pointsInARA.rawSparseData.ARAindex(jj)) %Then it's white matter
                    numWhiteMatter(ii)=numWhiteMatter(ii)+1;
                    if regexp(nType,'callosal','once')
                        numCallosal(ii)=numCallosal(ii)+1;
                        if ~any(wInds==data(ii).pointsInARA.rawSparseData.ARAindex(jj))
                            fprintf('** Neuron %d has a connection labeled as callosal which is not in the white matter\n',ii)
                        end
                        if ~isempty(regexp(nType,'bright','once'))
                            %bright callosal
                            numBrightCallosal(ii)=numBrightCallosal(ii)+1;
                        end
                    else
                        %Then it's white matter but not callosal
                        numWhiteMatterNonCallosal(ii)=numWhiteMatterNonCallosal(ii)+1;
                    end
                else
                    %It's gray matter
                    if ~isempty(regexp(nType,'fading','once'))
                        numFadingGrayMater(ii)=numFadingGrayMater(ii)+1;
                    end
                end


            end %for jj=1:length(oTrace.Node)

            if numAbruptTerminations(ii)==0 || ...
                (numAbruptTerminations(ii) == numCallosal(ii)) %So cells with failures only in the callosum get through
                cleanCell(ii)=1;
            end
            allCells(ii)=1;
        end %for ii=1:length(data)


        out.numberOfNeurons=n;
        out.V1only=V1only;
        out.numAbruptTerminations=sum(numAbruptTerminations);
        out.cellsWithAtLeastOneAbruptTermination = sum(numAbruptTerminations>0);
        out.numCallosal=sum(numCallosal>0);
        out.numBright=sum(numBright>0);
        out.numBrightCallosal=sum(numBrightCallosal>0);
        out.numFading=sum(numFading>0);
        out.numWhiteMatter=sum(numWhiteMatter>0);
        out.numWhiteMatterNonCallosal=sum(numWhiteMatterNonCallosal>0);
        out.numFadingGrayMater=sum(numFadingGrayMater);
        out.cleanCell=cleanCell;
        out.allCells=allCells;

        d=[data.details];

        out.indexValuesOfcleanCells = find(cleanCell);
        out.IDsofCleanCells = {d(find(cleanCell)).cellID};
        out.indexValuesOfAllCells = find(allCells);
        out.IDsOfAllCells = {d(find(allCells)).cellID};

function prettyDisplay(countData)
    if countData.V1only==1
        fprintf('\nV1 cells only -- ')
    elseif countData.V1only==-1
        fprintf('\nnon-V1 cells only -- ')
    elseif countData.V1only==0
        fprintf('\nBoth V1 and non V1 cells -- ')
    elseif isempty(countData.V1only)
        return
    end

    fprintf('%s\n',countData.backlabelled)

    fprintf('Number of cells: %d (%d have no abrupt terminations that aren''t callosal)\n', countData.numberOfNeurons, sum(countData.cleanCell));
    fprintf('Total abrupt terminations: %d\n', countData.numAbruptTerminations);
    fprintf('Cells with at least one abrupt termination: %d\n', countData.cellsWithAtLeastOneAbruptTermination);
    fprintf('Number of bright terminations: %d\n', countData.numBright);
    fprintf('Number of fading terminations: %d\n', countData.numFading);
    fprintf('Number of fading terminations in gray matter: %d\n', countData.numFadingGrayMater);
    fprintf('Number of abrupt terminations in white matter: %d\n', countData.numWhiteMatter);
    fprintf('Number of abrupt terminations in callosum: %d\n', countData.numCallosal);
    fprintf('Number of abrupt terminations in callosum that are bright: %d\n', countData.numBrightCallosal);
    fprintf('Number of abrupt terminations in white matter but not callosum: %d\n', countData.numWhiteMatterNonCallosal);


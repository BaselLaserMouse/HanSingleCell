function upSampledPoints = upsampleNeuriteTree(neuriteTree,resampleAtInterval)
% Upsample a neurite tree and return all points as a matrix (so not as a tree)  
%
% 
%   function upSampledPoints = upsampleNeuriteTree(neuriteTree,resampleAtInterval)
%
% 
% Purpose
% Upsample a traced neurite tree such that we have one sample per unit. 
% e.g. if the tree was downsampled for 25 micron voxels (or whatever voxel size our
% atlas had) and we set resampleAtInterval to 1, then we get one point per voxel. So one 
% every 25 microns, if we're working with the 25 micron tree. 
% If we set resample to, say 0.5 and we have 25 micron voxels then we get one point every 25*0.5 = 12.25 microns
%
%
% Inputs
% neuriteTree - this is the output of exportedCSV2tree
% resampleAtInterval - number of points per voxel
% 
%
% Outputs
% upsampledPoints - The output is an array where each row is one index and the 
%   columns are: z, x, and y. This can be read by pointsInARA. We don't return 
%   the complete upsampled tree, as there's no point. We're just going to use 
%   these data to estimate neurite length by brain area.
%
% 
%
% Rob Campbell - Basel 2016
%
%
% Also see: aratools.sparse.processTreeData



segments = neuriteTree.getsegments;


upSampledPoints=[];

for ii=1:length(segments)
    thisSegment = segments{ii};
    if length(thisSegment)<2
        continue
    end

    pointsFromSegment = [neuriteTree.Node{thisSegment}];
    pointsFromSegment = reshape(pointsFromSegment, 3, [])';

    upSampledPoints = [upSampledPoints; doUpSampleOfSegment(pointsFromSegment,resampleAtInterval)];
end




% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
function upsampled = doUpSampleOfSegment(points,resampleAtInterval)
    %upsample a segment
    % 
    % Inputs
    % points - n by 3 matrix of points from one axon segment
    % resampleAtinterval - number of points per voxel to produce
    %
    % Outputs
    % upampled - the returned upsampled tree

    Z=[];
    X=[];
    Y=[];

    for ii=1:size(points,1)-1
        ind = ii:ii+1;

        %Go through the data set point by point to interpolate between each pair of points. 
        [z,x,y] = doResample(points(ind,1), points(ind,2), points(ind,3), resampleAtInterval);

        Z=[Z;z];
        X=[X;x];
        Y=[Y;y];
    end

    upsampled = [Z,X,Y];




% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
function [rZ,rX,rY] = doResample(Z,X,Y, resampleAtInterval)
    %Upsample between two points defined in 3D space (The three inputs are vectors of length 2)
    d = sqrt( diff(Z)^2 + diff(X)^2 + diff(Y)^2);

    if d>resampleAtInterval
        %we will resample into this many equally spaced steps
        n=round(d/resampleAtInterval);

        rZ = linspace(Z(1),Z(2),n)';
        rX = linspace(X(1),X(2),n)';
        rY = linspace(Y(1),Y(2),n)';
    
        %check
        d = sqrt( diff(rZ).^2 + diff(rX).^2 + diff(rY).^2);
        if std(d)>resampleAtInterval*0.1
            fprintf('\n\nTHERE MAY BE AN ERROR IN THE RESAMPLING! -- The distance between points is variable\n\n')
        end
    else %don't attempt to resample if we're already smaller than the distance between points
        fprintf('Skipping segment resample. Points are already finer resolution in this segment\n')
        rZ=Z;
        rX=X;
        rY=Y;
    end




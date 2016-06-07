function [skelData, varargout] = loadSkeleta(filename,startIndcs,stopIndcs,minLength,verbose)
% go through data set, load data in chunks, and assemble full skeleton
% array from that
if nargin<4
    minLength = 1;
end
if nargin<5
    verbose = false;
end
% read out the size of the dataset
[~,h5Dims] = H5S.get_simple_extent_dims(H5D.get_space(H5D.open(...
    H5F.open(filename),'/skeleton')));
% NOTE: the long dimension is the first one of h5dims, but the last one of
% the dataset

regionSizes = stopIndcs - startIndcs + 1;
% only keep regions that are bigger than at least minLength
if minLength>1
    keepIdcs = regionSizes>=minLength;
    regionSizes = regionSizes(keepIdcs);
    startIndcs = startIndcs(keepIdcs);
    stopIndcs = stopIndcs(keepIdcs);
end

numRegions = length(startIndcs);
numData = sum(regionSizes);
skelData = NaN(2,49,numData);
skelDataIndcs = cumsum([1; regionSizes]);
if nargout > 1
    % load meta data
    trajectoryData = h5read(filename,'/trajectories_data');
    wormIDs = NaN(numData,1);
    if nargout > 2
        frameIDS = NaN(numData,1);
    end
end

for regCtr = 1:numRegions
    % set frame up to which to read, or end of data set
    upToFrame = stopIndcs(regCtr);
    fromFrame = startIndcs(regCtr);
    % load skeletal data
    if verbose
        display(['reading frames ' num2str(fromFrame) ' to ' num2str(upToFrame)...
            ' out of ' num2str(h5Dims(1)) ' frames for ' filename(end-42:end-5)])
    end
    skelData(:,:,skelDataIndcs(regCtr):(skelDataIndcs(regCtr + 1) - 1)) =...
        h5read(filename,'/skeleton',[1 1 fromFrame],...
        [2, 49, regionSizes(regCtr)]);
    if nargout > 1
        wormIDs(skelDataIndcs(regCtr):(skelDataIndcs(regCtr + 1) - 1)) = ...
            trajectoryData.worm_index_joined(fromFrame:upToFrame);
        if nargout > 2
            frameIDs(skelDataIndcs(regCtr):(skelDataIndcs(regCtr + 1) - 1)) = ...
                trajectoryData.frame_number(fromFrame:upToFrame);
        end
    end
end

if nargout > 1
    varargout{1} = wormIDs;
    if nargout > 2
       varargout{2} = frameIDs; 
    end
end
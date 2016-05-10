function skelData = loadSkeleta(filename,startIndcs,stopIndcs,verbose)
% go through data set, load data in chunks, and assemble full skeleton
% array from that

% read out the size of the dataset
[~,h5Dims] = H5S.get_simple_extent_dims(H5D.get_space(H5D.open(...
    H5F.open(filename),'/skeleton')));
% NOTE: the long dimension is the first one of h5dims, but the last one of
% the dataset

regionSizes = stopIndcs - startIndcs + 1;
numRegions = length(startIndcs);
skelData = NaN(2,49,sum(regionSizes));
skelDataIndcs = cumsum([1; regionSizes]);

for regCtr = 1:numRegions
    % set frame up to which to read, or end of data set
    upToFrame = stopIndcs(regCtr);
    fromFrame = startIndcs(regCtr);
    % load skeletal data
    if verbose
        display(['reading frames ' num2str(fromFrame) ' to ' num2str(upToFrame)...
            ' out of ' num2str(h5Dims(1)) ' frames for ' filename(48:end-5)])
    end
    skelData(:,:,skelDataIndcs(regCtr):(skelDataIndcs(regCtr + 1) - 1)) =...
        h5read(filename,'/skeleton',[1 1 fromFrame],...
        [2, 49, regionSizes(regCtr)]);   
end
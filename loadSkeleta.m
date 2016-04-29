function skelData = loadSkeleta(filename,nFrames,verbose)
% go through data set, load data in chunks, and assemble full skeleton
% array from that

% read out the size of the dataset
[~,h5Dims] = H5S.get_simple_extent_dims(H5D.get_space(H5D.open(...
    H5F.open(filename),'/skeleton')));
% NOTE: the long dimension is the first one of h5dims, but the last one of
% the dataset

skelData = NaN(2,49,h5Dims(1));

for frameCtr = 1:nFrames:h5Dims(1)
    % set frame up to which to read, or end of data set
    upToFrame = min(frameCtr + nFrames - 1,h5Dims(1));
    % load skeletal data
    if verbose
        display(['reading frames ' num2str(frameCtr) ' to ' num2str(upToFrame)...
            ' out of ' num2str(h5Dims(1)) ' frames for ' filename(48:end-5)])
    end
    
    skelData(:,:,frameCtr:upToFrame) = h5read(filename,'/skeleton',[1 1 frameCtr],...
        [2, 49, (upToFrame - frameCtr + 1)]);   
end
function out = filterData(filename,nFrames,verbose)
% filter worm skeletal data by area, observed time, etc.
% speed filtering is hierarchical, meaning that the first few filters will 
% reduce the number of data to pass through the speed filter, to reduce 
% computational intensity of computing speeds for every worm

% load all metadata
trajectoryData = h5read(filename,'/trajectories_data');

% select frames that have skeleton
hasSkel = trajectoryData.has_skeleton==1;

% select frames with worms that occurr more than a certain number
framesPerWorm = histcounts(trajectoryData.worm_index_joined,max(trajectoryData.worm_index_joined));
frequentWorms = find(framesPerWorm>=25*10);
framesFilter = ismember(trajectoryData.worm_index_joined,frequentWorms);

% select frames with a certain area
areaFilter = filterArea(trajectoryData,25,1500,50,50,hasSkel&framesFilter);

% select worms with at least a certain speed
speedFilter = filterSpeed(trajectoryData,0.1,hasSkel&framesFilter&areaFilter);

% detect dust from labelled data
dustIdcs = filterDust(trajectoryData,0.1,5);

% load skeletal data
skelData = loadSkeleta(filename,nFrames,verbose);

% filter data
out = skelData(:,:,hasSkel&areaFilter&framesFilter&speedFilter&~dustIdcs);
end
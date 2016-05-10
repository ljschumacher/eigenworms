function out = filterData(filename,verbose)
% filter worm skeletal data by area, observed time, etc.
% some filtering is hierarchical, meaning that the first few filters will 
% reduce the number of data to pass through the next filter

% issues/to-do:
% - should we only load data patches of a minimum contiguous length? (after
% having already filtered for how many frames a worm appears in, which can
% be reduced by other filters)

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

% detect dust from manually labelled data
dustIdcs = filterDust(trajectoryData,0.1,5);

% combine filters and select data regions to load
combiFilter = hasSkel&areaFilter&framesFilter&speedFilter&~dustIdcs;
dF = diff([0; combiFilter; 0]);
startIndcs = find(dF==1); % pick out start of contiguous data regions
stopIndcs = find(dF==-1) - 1; % pick out stop of contiguous data regions

% load skeletal data
skelData = loadSkeleta(filename,startIndcs,stopIndcs,verbose);

% filter data
out = skelData;
end
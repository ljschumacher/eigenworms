function [skelData, varargout] = filterData(filename,verbose,plotDiagnostics,minFrames)
% filter worm skeletal data by area, observed time, etc.
% some filtering is hierarchical, meaning that the first few filters will
% reduce the number of data to pass through the next filter

if nargin <4
    minFrames = 25*30;
end

if verbose
    display(['now filtering dataset ' filename(end-42:end-5)])
end

% load all metadata
trajectoryData = h5read(filename,'/trajectories_data');

% select frames that have skeleton
hasSkel = trajectoryData.has_skeleton==1;

% select frames with worms that occurr more than a certain number
framesPerWorm = histcounts(trajectoryData.worm_index_joined,max(trajectoryData.worm_index_joined));
frequentWorms = find(framesPerWorm>=minFrames);
frequentFilter = ismember(trajectoryData.worm_index_joined,frequentWorms);

% select frames with a certain area
areaFilter = filterArea(trajectoryData,20,50,hasSkel&frequentFilter,...
    plotDiagnostics,filename);

% detect dust from manually labelled data
dustIdcs = filterDust(trajectoryData,0.1,5);

% select worms with at least a certain speed
speedFilter = filterSpeed(trajectoryData,0.1,verbose,...
    hasSkel&frequentFilter&areaFilter&~dustIdcs,plotDiagnostics,filename);

% combine filters and select data regions to load
combiFilter = hasSkel&areaFilter&frequentFilter&speedFilter&~dustIdcs;
dFrame = diff([0; combiFilter; 0]);
startIndcs = find(dFrame==1); % pick out start of contiguous data regions
stopIndcs = find(dFrame==-1) - 1; % pick out stop of contiguous data regions

if plotDiagnostics
    plotWormNumbers(filename,trajectoryData,combiFilter)
end
% load  (filtered) skeletal data
switch nargout
    case 1
        skelData = loadSkeleta(filename,startIndcs,stopIndcs,minFrames,verbose);
    case 2
        [skelData, varargout{1}] = loadSkeleta(filename,startIndcs,stopIndcs,minFrames,verbose);
    case 3
        [skelData, varargout{1}, varargout{2}] = loadSkeleta(filename,startIndcs,stopIndcs,minFrames,verbose);
end
end
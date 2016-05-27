function movingIdcs = filterSpeed(trajectoryData,speedThreshold,verbose,otherFilters)
% find all worms that move at least once above a given threshold speed
if nargin<4
    otherFilters = true(size(trajectoryData.has_skeleton));
end
wormIDs = unique(trajectoryData.worm_index_joined(otherFilters))';
movingIdcs = false(size(trajectoryData.frame_number)); % initialise all frames as excluded
for wormCtr=1:length(wormIDs) % go through worms and calculate speeds
    if verbose
        display(['calculating speed for tracked object ' num2str(wormCtr) ...
            ' out of ' num2str(length(wormIDs))])
    end
    % find all frames the current worm is in, excluding previously filtered
    % frames to improve performance (hopefully)
    wormIdcs = trajectoryData.worm_index_joined==wormIDs(wormCtr)&otherFilters;
    wormDx = diff(trajectoryData.coord_x(wormIdcs));
    wormDy = diff(trajectoryData.coord_y(wormIdcs));
    wormDf = diff(trajectoryData.frame_number(wormIdcs));
    wormDisplacement = sqrt(wormDx.^2 + wormDy.^2)./single(wormDf);
    movingIdcs(wormDisplacement>speedThreshold) = 1; % include worm-frames were movement is above threshold
end
function movingIdcs = filterSpeed(trajectoryData,speedThreshold)
% find all worms that move at least once above a given threshold speed
wormIDs = unique(trajectoryData.worm_index_joined)';
wormMoving = false(size(wormIDs));
for wormCtr=1:length(wormIDs) % go through worms and calculate speeds
    wormIdcs = trajectoryData.worm_index_joined==wormIDs(wormCtr);
    wormDx = diff(trajectoryData.coord_x(wormIdcs));
    wormDy = diff(trajectoryData.coord_y(wormIdcs));
    wormDf = diff(trajectoryData.frame_number(wormIdcs));
    wormDisplacement = sqrt(wormDx.^2 + wormDy.^2)./single(wormDf);
    wormMoving(wormCtr) = any(wormDisplacement>speedThreshold);
end
movingWorms = wormIDs(wormMoving);
movingIdcs = ismember(trajectoryData.worm_index_joined,movingWorms);
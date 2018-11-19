function dustIdcs = filterDust(trajectoryData,relativeAreaTolerance,relativePositionTolerance)
% detect dust from 'bad' worm labels and find objects with similar area and
% position
if isfield(trajectoryData,'worm_label') % filter for worm labels, if they exist
    % find the frames with bad/dust stuffs
    dustIdcs = trajectoryData.worm_label==3;
    dustIDs=unique(trajectoryData.worm_index_joined(dustIdcs));
    % also remove non-labelled dust elements (according to area and x and y positions)
    for wormCtr=1:numel(dustIDs)
        meanArea = mean(trajectoryData.area(dustIdcs));
        relativeDiffArea=abs(trajectoryData.area - meanArea)/meanArea;
        relativeDiffX=abs(trajectoryData.coord_x - mean(trajectoryData.coord_x(dustIdcs)));
        relativeDiffY=abs(trajectoryData.coord_y - mean(trajectoryData.coord_y(dustIdcs)));
        
        dustIdcs(relativeDiffArea<relativeAreaTolerance&...
            relativeDiffX<relativePositionTolerance&...
            relativeDiffY<relativePositionTolerance) = true;
    end
else
    display(['No worm_label found for this dataset'])
    dustIdcs = false(size(trajectoryData.worm_index_joined));
end

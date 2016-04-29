function areaIdcs = filterArea(trajectoryData,binWidth,binRange,...
    minPeakWidth,minPeakDistance,otherFilters)
% identify the most prominent peak in the area distribution
bins = 0:binWidth:binRange;
counts = histcounts(trajectoryData.area(otherFilters),bins,'normalization','pdf');
[~, locs, widths, proms] = findpeaks(counts,bins(2:end)-binWidth/2,...
    'MinPeakDistance',minPeakDistance,'MinPeakWidth',minPeakWidth);
[~, mostProm] = max(proms); % find most prominent peak
areaIdcs = trajectoryData.area>=(locs(mostProm) - widths(mostProm))&...
    trajectoryData.area<=(locs(mostProm) + widths(mostProm));
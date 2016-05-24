function areaIdcs = filterArea(trajectoryData,binWidth,binRange,...
    minPeakWidth,minPeakDistance,otherFilters,plotDiagnostics,filename)
% identify the most prominent peak in the area distribution
if nargin<6
    otherFilters = true(size(trajectoryData.has_skeleton));
end
if nargin < 7
    plotDiagnostics = 0;
    filename = [];
end

bins = 0:binWidth:binRange;

if plotDiagnostics
    areaHistFig = figure;
    h = histogram(trajectoryData.area(otherFilters),bins,'normalization','pdf');
    hold on
    [peaks, locs, widths, proms] = findpeaks(h.Values,bins(2:end)-binWidth/2,...
        'MinPeakDistance',minPeakDistance,'MinPeakWidth',minPeakWidth);
    plot(locs,peaks,'v','LineWidth',2)
    [~, mostProm] = max(proms); % find most prominent peak
    plot(locs(mostProm),peaks(mostProm),'ro','LineWidth',2,'MarkerSize',20)
    stairs([locs - widths; locs - widths; locs + widths],...
        [0; 1; 0]*peaks,'LineWidth',2)
    xlim([0 2500])
    xlabel('area of tracked opbject')
    ylabel('pdf')
    % save plot
    figName = ['figures/diagnostics/areaHist_dataset_' filename(end-42:end-15) '.eps'];
    exportfig(gcf,figName,'Color','rgb')
    system(['epstopdf ' figName]);
    % clean up
    system(['rm ' figName]);
    close(areaHistFig)
else
    counts = histcounts(trajectoryData.area(otherFilters),bins,'normalization','pdf');
    [~, locs, widths, proms] = findpeaks(counts,bins(2:end)-binWidth/2,...
        'MinPeakDistance',minPeakDistance,'MinPeakWidth',minPeakWidth);
    [~, mostProm] = max(proms); % find most prominent peak
end
areaIdcs = trajectoryData.area>=(locs(mostProm) - widths(mostProm))&...
    trajectoryData.area<=(locs(mostProm) + widths(mostProm));
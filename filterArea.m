function areaIdcs = filterArea(trajectoryData,...
    minPeakWidth,minPeakDistance,otherFilters,plotDiagnostics,filename)
% identify the most prominent peak in the area distribution
if nargin < 4
    otherFilters = true(size(trajectoryData.has_skeleton));
end
if nargin < 5
    plotDiagnostics = 0;
    filename = [];
end

if plotDiagnostics
    areaHistFig = figure;
    h = histogram(trajectoryData.area(otherFilters),'BinMethod','sqrt','normalization','pdf',...
        'EdgeColor','none');
    hold on
    [peaks, locs, widths, proms] = findpeaks(h.Values,double(h.BinEdges(2:end)-h.BinWidth/2),...
        'MinPeakDistance',minPeakDistance,'MinPeakWidth',max(minPeakWidth,2*h.BinWidth));
    plot(locs,peaks,'v','LineWidth',2)
    [~, mostProm] = max(proms); % find most prominent peak
    plot(locs(mostProm),peaks(mostProm),'ro','LineWidth',2,'MarkerSize',20)
    stairs([locs - widths; locs - widths; locs + widths],...
        [0; 1; 0]*peaks,'LineWidth',2)
    xlabel('area of tracked object')
    ylabel('pdf')
    title(filename(end-42:end-15),'Interpreter','none')
    % save plot
    figName = ['figures/diagnostics/areaHist_dataset_' filename(end-42:end-15) '.eps'];
    exportfig(gcf,figName,'Color','rgb','LockAxes',0); % lock axes 0 is important for keeping 10^x axes labels
    system(['epstopdf ' figName]);
    % clean up
    system(['rm ' figName]);
    close(areaHistFig)
else
    [counts, bins] = histcounts(trajectoryData.area(otherFilters),'normalization','pdf');
    binWidth = mean(diff(bins));
    [~, locs, widths, proms] = findpeaks(counts,double(bins(2:end)-binWidth/2),...
        'MinPeakDistance',minPeakDistance,'MinPeakWidth',max(minPeakWidth,2*binWidth));
    [~, mostProm] = max(proms); % find most prominent peak
end
areaIdcs = trajectoryData.area>=(locs(mostProm) - widths(mostProm))&...
    trajectoryData.area<=(locs(mostProm) + widths(mostProm));
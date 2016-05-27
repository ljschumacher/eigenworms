function movingIdcs = filterSpeed(trajectoryData,speedThreshold,...
    verbose,otherFilters,plotDiagnostics,filename)
% find all worms that move at least once above a given threshold speed
if nargin<4
    otherFilters = true(size(trajectoryData.has_skeleton));
end
if nargin<5
    plotDiagnostics = 0;
    filename = [];
end
wormIDs = unique(trajectoryData.worm_index_joined(otherFilters))';
nWorms = length(wormIDs);
movingIdcs = false(size(trajectoryData.frame_number)); % initialise all frames as excluded
if plotDiagnostics
    dBin = 0.05;
    bins = 0:dBin:2.5;
    speedDistributions = NaN(nWorms,length(bins)-1);
end
for wormCtr=1:nWorms % go through worms and calculate speeds
    if verbose
        display(['calculating speed for tracked object ' num2str(wormCtr) ...
            ' out of ' num2str(nWorms)])
    end
    % find all frames the current worm is in, excluding previously filtered
    % frames to improve performance (hopefully)
    wormIdcs = trajectoryData.worm_index_joined==wormIDs(wormCtr)&otherFilters;
    wormDx = diff(trajectoryData.coord_x(wormIdcs));
    wormDy = diff(trajectoryData.coord_y(wormIdcs));
    wormDf = diff(trajectoryData.frame_number(wormIdcs));
    wormSpeed = sqrt(wormDx.^2 + wormDy.^2)./single(wormDf);
    movingIdcs(wormSpeed>speedThreshold) = 1; % include worm-frames were movement is above threshold
    if plotDiagnostics
        speedDistributions(wormCtr,:) = histcounts(wormSpeed,bins,'Normalization','pdf');
    end
end
if plotDiagnostics
    speedHistFig = figure;
    plot(bins(1:end-1)+dBin/2,speedDistributions')
    hold on
    plot([speedThreshold, speedThreshold], [0, max(max(speedDistributions))],...
        'k--','LineWidth',2)
    title(filename(end-42:end-15),'Interpreter','none')
    ylabel('pdf'), xlabel('speed (pixel/frame)')
    % save plot
    figName = ['figures/diagnostics/speedHist_dataset_' filename(end-42:end-15) '.eps'];
    exportfig(gcf,figName,'Color','rgb');
    system(['epstopdf ' figName]);
    % clean up
    system(['rm ' figName]);
    close(speedHistFig)
end
function movingIdcs = filterSpeed(trajectoryData,speedThreshold,...
    verbose,otherFilters,plotDiagnostics,filename)
% find all worms that move at least once above a given threshold speed
if nargin<3
    verbose = 0; % quiet by default
end
if nargin<4
    otherFilters = true(size(trajectoryData.has_skeleton)); % no other filtering by default
end
if nargin<5
    plotDiagnostics = 0; % no diagnostic plots by default
    filename = [];
end
wormIDs = unique(trajectoryData.worm_index_joined(otherFilters))';
nWorms = length(wormIDs);
movingIdcs = false(size(trajectoryData.frame_number)); % initialise all frames as excluded
if verbose, displayEvery = 100; end % set how often to display status
if plotDiagnostics % set parameters for calculating speed distributions
    dBin = 0.05;
    bins = 0:dBin:2.5;
    speedDistributions = NaN(nWorms,length(bins)-1);
end
for wormCtr=1:nWorms % go through worms and calculate speeds
    if verbose&&~mod(wormCtr,displayEvery)
        display(['calculating speed for tracked objects ' num2str(wormCtr) ...
            ' to ' num2str(wormCtr+displayEvery) ' out of ' num2str(nWorms)])
    end
    % find all frames the current worm is in, excluding previously filtered
    % frames to improve performance (hopefully)
    wormIdcs = trajectoryData.worm_index_joined==wormIDs(wormCtr)&otherFilters;
    wormDx = diff(trajectoryData.coord_x(wormIdcs));
    wormDy = diff(trajectoryData.coord_y(wormIdcs));
    wormDf = diff(trajectoryData.frame_number(wormIdcs));
    wormSpeed = sqrt(wormDx.^2 + wormDy.^2)./single(wormDf);
    frameNums = find(wormIdcs); % refer back to frame numbers in dataset
    movingIdcs(frameNums(wormSpeed>speedThreshold)) = 1; % include worm-frames were movement is above threshold
    if plotDiagnostics % calculate distribution of speeds for this worm
        speedDistributions(wormCtr,:) = histcounts(wormSpeed,bins,'Normalization','pdf');
    end
end
if plotDiagnostics % plot speed distribitions
    speedHistFig = figure;
    % fourth element of color vector sets transparancy
    plot(bins(1:end-1)+dBin/2,speedDistributions','Color',[0 0 1 1/nWorms])
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
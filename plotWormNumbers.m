function plotWormNumbers(filename,trajectoryData,filter)

% figure export options
exportOptions = struct('Color','rgb','LockAxes',0); % lock axes 0 is important for keeping 10^x axes labels

dFrame = 25*4; % only plot every nFrame-th frame

% plot worm numbers pre-filtering
frameNums = min(trajectoryData.frame_number):dFrame:max(trajectoryData.frame_number);
framesFilter = ismember(trajectoryData.frame_number,frameNums);
wormNumFig = figure;
histogram(trajectoryData.frame_number(framesFilter),frameNums,'DisplayStyle','stairs');
hold on
% combine filters
combiFilter = filter&framesFilter;
if nnz(combiFilter)
    % plot post-filtered worm numbers
    filteredFrameNums = unique(trajectoryData.frame_number(combiFilter))';
    [wormNums, wormBins] = histcounts(trajectoryData.frame_number(combiFilter),filteredFrameNums);
    stairs(wormBins(1:end-1),wormNums); % don't use histogram again to make it parfor compatible
else
    display('We filtered out all the worms! :(')
end
xlabel('frame number')
ylabel('object count')
legend('raw','filtered')
xlim([0 max(frameNums)])
title(filename(end-42:end-15),'Interpreter','none')
% save plot
figName = ['figures/diagnostics/wormNums_dataset_' filename(end-42:end-15) '.eps'];
exportfig(gcf,figName,exportOptions)
system(['epstopdf ' figName]);
% clean up
system(['rm ' figName]);
close(wormNumFig)
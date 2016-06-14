function plotShapeCorrelations(timeSeries,wormIDarray,...
    frameIDarray,maxLag,savepath,componentWise)
% calculate and plot auto correlations of projected amplitudes in worm data
% INPUTS
% timeSeries - this should be an angle array or projected amplitudes of
% shape principal components, shape is m-wormframes by n-dimensions
% wormIDarray - vector of same length as timeSeries with a unique number
% identifiying the track
% frameIDarray - optional (default assumes time-ordered tracks) vector of
% same length as timeSeries with a unique number for to identify temporal
% order within a track
% maxLag - optional (default 25) - maximum lag (in number of frames) for
% autocorrelation calculation
% savepath - optional (default current directory), path where figure is to
% be saved
% componentWise - optional (default false), specifies whether to compute
% autocorrelation for each component or via the scalar product

% issues/to-do:
% - should we normalise the variance of different tracks?
% - spacing of time series may not be uniform, so lag could be inconcistent

% figure export options
exportOptions = struct('Color','rgb');

if nargin<3
    needSort = 'False';
    if nargin<4
        maxLag = 25;
        if nargin<5
            savepath = './';
            if nargin<6
                componentWise = 0;
            end
        end
    end
else
    needSort = 'True';
end

pacFig = figure;

% loop through worms
wormIDs = unique(wormIDarray);
nWorms = length(wormIDs);
nDims = size(timeSeries,2);
if componentWise
    projAutoCorr = NaN(nWorms,nDims,maxLag + 1);
else
    projAutoCorr = NaN(nWorms,maxLag + 1);
end
for wormCtr = 1:nWorms
    % find frames with current worm
    wormIdcs = wormIDarray==wormIDs(wormCtr);
    if needSort
        % select (unsorted) timeseries of current worm
        currentProjections = timeSeries(wormIdcs,:);
        % bring randomly sampled frame back into temporal order
        [~, sortInd] = sort(frameIDarray(wormIdcs));
        % calculate autocorrelation of (sorted) projections
        if componentWise
            for dimCtr = 1:nDims
                projAutoCorr(wormCtr,dimCtr,:) = vectorAutoCorrelation(currentProjections(sortInd,dimCtr)',maxLag);
            end
        else
            projAutoCorr(wormCtr,:) = vectorAutoCorrelation(currentProjections(sortInd,:)',maxLag);
        end
    else
        if componentWise
            for dimCtr = 1:nDims
                projAutoCorr(wormCtr,dimCtr,:) = vectorAutoCorrelation(timeSeries(wormIdcs,dimCtr)',maxLag);
            end
        else
            projAutoCorr(wormCtr,:) = vectorAutoCorrelation(timeSeries(wormIdcs,:)',maxLag);
        end
    end
end
if componentWise
    plot(squeeze(nanmean(projAutoCorr))')
    legend(num2str((1:nDims)'))
    title(['mean of ' num2str(nWorms) ' tracks, '...
        num2str(size(timeSeries,1)/25/3600,2) ' worm-hours, first '...
        num2str(nDims) ' shape dimensions'],'FontWeight','normal');
else
    plot(nanmean(projAutoCorr)) % nanmean may be important due to how the autocorr is calculated eg for long lags with fewer samples
    title(['mean of ' num2str(nWorms) ' tracks, '...
        num2str(size(timeSeries,1)/25/3600,2) ' worm-hours, '...
        num2str(nDims) ' shape dimensions'],'FontWeight','normal');
end
% annotate and save figure
ylabel('autocorrelation')
xlabel('lag (frames)')
xlim([0 maxLag])
figName = [savepath 'shape_autocorr.eps'];
exportfig(pacFig,figName,exportOptions)
system(['epstopdf ' figName]);
system(['rm ' figName]);
close(pacFig)
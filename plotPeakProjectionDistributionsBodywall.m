% plots the distribution of local extrema of eigenworm projections
close all
clear

% issues/to-do:
% - filter red channel data further? currently we are assuming that those
% with feature calculations are filtered fairly well
% - can we vectorize the loop over frames?

% figure export options
exportOptions = struct('Color','rg');

numEigenworms = 4; % number of projections to plot

neighbourCutOff = 500; % distance in microns to consider a neighbour close
minNumNeighbours = 2;
strains = {'npr1','N2'};
wormnums = {'40','HD','1W'};
intensityThresholds_g = [60, 40, NaN];

pixelsize = 100/19.5; % 100 microns are 19.5 pixels
maxBlobSize_g = 1e4; % for filtering green channel worms

for numCtr = 1:length(wormnums)
    wormnum = wormnums{numCtr};
    for strainCtr = 1:length(strains)
        strain = strains{strainCtr};
        %% load data
        filenames = importdata([strains{strainCtr} '_' wormnum '_r_list.txt']);
        numFiles = length(filenames);
        if ~strcmp(wormnum,'1W')
            filenames_g = importdata([strains{strainCtr} '_' wormnum '_g_list.txt']);
        else
            filenames_g = {};
        end
        peakProjections1 = cell(numFiles,1);
        peakProjections2 = cell(numFiles,1);
        peakProjections3 = cell(numFiles,1);
        peakProjections4 = cell(numFiles,1);
        for fileCtr = 1:numFiles
            if ~strcmp(wormnum,'1W'), filename_g = filenames_g{fileCtr}; end
            filename = filenames{fileCtr};
            if exist(strrep(filename,'skeletons','features'),'file')&&exist(filename,'file')...
                    &&(strcmp(wormnum,'1W')||exist(filename_g,'file'))
                features = h5read(strrep(filename,'skeletons','features'),'/features_timeseries');
                trajData = h5read(filename,'/trajectories_data');
                if ~strcmp(wormnum,'1W')
                    trajData_g = h5read(filename_g,'/trajectories_data');
                    blobFeats_g = h5read(filename_g,'/blob_features');
                end
                if ~strcmp(wormnum,'1W') % if it is multiworm data, we need to filter for worms in clusters
                    features.filtered = false(size(features.timestamp));
                    framesAnalyzed = unique(features.timestamp);
                    numFrames = numel(framesAnalyzed);
                    % filter green channel by blob size and intensity
                    trajData_g.filtered = (blobFeats_g.area*pixelsize^2<=maxBlobSize_g)&...
                        (blobFeats_g.intensity_mean>=intensityThresholds_g(numCtr));
                    % calculate red-green neighbour distances to filter for in-cluster
                    for frameCtr = 1:numFrames
                        frame = framesAnalyzed(frameCtr);
                        % as there are more objects in trajData than
                        % features, we need to select the right subset
                        trajData.filtered = ismember(trajData.worm_index_joined,...
                            int32(features.worm_index(features.timestamp==frame)));
                        [x, y] = getWormPositions(trajData, frame);
                        [x_g, y_g] = getWormPositions(trajData_g, frame);
                        if numel(x_g)>=1&&numel(x)>=1 % need at least two worms in frame to calculate distances
                            redToGreenDistances = pdist2([x y],[x_g y_g]).*pixelsize; % distance of every red worm to every green
                            numNeighbours = sum(redToGreenDistances<neighbourCutOff,2);
                        elseif numel(x)>=1
                            numNeighbours = zeros(size(x));
                        else
                            numNeighbours = 0; % sets all frames to false if there are no valid worms in current frame
                        end
                        % exclude worm-frames with fewer than some number of neighbours
                        features.filtered(features.timestamp==frame) = ...
                            numNeighbours>=minNumNeighbours; % keeps only those worms more than minNumNeighbours neighbours
                    end
                else
                    features.filtered = true(size(features.timestamp));
                end
                %% find peaks and troughs in eigenworm projections
                % exclude when the worm index changes, to break continuity of time series
                features.filtered(diff(features.worm_index)~=0) = false;
                if any(features.filtered)
                    features.eigen_projection_1 = double(features.eigen_projection_1);
                    features.eigen_projection_2 = double(features.eigen_projection_2);
                    features.eigen_projection_3 = double(features.eigen_projection_3);
                    features.eigen_projection_4 = double(features.eigen_projection_4);
                    peakProjections1{fileCtr} = [findpeaks(...
                        features.eigen_projection_1(features.filtered),'MinPeakProminence',0.5);
                        -findpeaks(-features.eigen_projection_1(features.filtered),'MinPeakProminence',0.5)];
                    peakProjections2{fileCtr} = [findpeaks(...
                        features.eigen_projection_2(features.filtered),'MinPeakProminence',0.5);
                        -findpeaks(-features.eigen_projection_2(features.filtered),'MinPeakProminence',0.5)];
                    peakProjections3{fileCtr} = [findpeaks(...
                        features.eigen_projection_3(features.filtered),'MinPeakProminence',0.5);
                        -findpeaks(-features.eigen_projection_3(features.filtered),'MinPeakProminence',0.5)];
                    peakProjections4{fileCtr} = [findpeaks(...
                        features.eigen_projection_4(features.filtered),'MinPeakProminence',0.5);
                        -findpeaks(-features.eigen_projection_4(features.filtered),'MinPeakProminence',0.5)];
                else
                    warning(['All worms filtered out for ' filename ])
                end
            else
                warning(['Not all necessary tracking results present for ' filename ])
            end
        end
        % pool data from multiple recordings
        peakProjections1 = vertcat(peakProjections1{:});
        peakProjections2 = vertcat(peakProjections2{:});
        peakProjections3 = vertcat(peakProjections3{:});
        peakProjections4 = vertcat(peakProjections4{:});
        %% plot data
        subplot(4,2,strainCtr)
        if numCtr==1
            hold on, title(strain,'FontWeight','normal')
            ylabel('P'), ylim([0 0.28]), yticks(get(gca,'ylim')), xlim([-10 10]), xlabel('a^*_1'), box on
        end
        histogram(peakProjections1,'Normalization','pdf','DisplayStyle','stairs')
 
        subplot(4,2,strainCtr+2)
        if numCtr==1
            hold on, ylabel('P'), ylim([0 0.28]), yticks(get(gca,'ylim')), xlim([-10 10]), xlabel('a^*_2'), box on
        end
        histogram(peakProjections2,'Normalization','pdf','DisplayStyle','stairs')
        
        subplot(4,2,strainCtr+4)
        if numCtr==1
            hold on, ylabel('P'), ylim([0 0.28]), yticks(get(gca,'ylim')), xlim([-10 10]), xlabel('a^*_3'), box on
        end
        histogram(peakProjections3,'Normalization','pdf','DisplayStyle','stairs')

        subplot(4,2,strainCtr+6)
        if numCtr==1
            hold on, ylabel('P'), ylim([0 0.28]), yticks(get(gca,'ylim')), xlim([-10 10]), xlabel('a^*_4'), box on
        end
        histogram(peakProjections4,'Normalization','pdf','DisplayStyle','stairs')
    end
end
%% format and export figure
set(gcf,'Position',get(gcf,'Position').*[1 1 1 1.5]) % make figure higher
legend(wormnums)
figFileName = ['figures/peakAmplitudeEigenwormProjectionDistribution.eps'];
exportfig(gcf,figFileName,'Color','rgb')
system(['epstopdf ' figFileName]);
system(['rm ' figFileName]);
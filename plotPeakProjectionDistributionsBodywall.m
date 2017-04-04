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

minNumNeighbours = 3;
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
        filenames = importdata(['datalists/' strains{strainCtr} '_' wormnum '_r_list.txt']);
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
            filename = filenames{fileCtr};
            if ~strcmp(wormnum,'1W'), filename_g = filenames_g{fileCtr}; end
            if exist(strrep(filename,'skeletons','features'),'file')&&exist(filename,'file')...
                    &&(strcmp(wormnum,'1W')||exist(filename_g,'file'))
                features = h5read(strrep(filename,'skeletons','features'),'/features_timeseries');
                trajData = h5read(filename,'/trajectories_data');
                if ~strcmp(wormnum,'1W')
                    trajData_g = h5read(filename_g,'/trajectories_data');
                    blobFeats_g = h5read(filename_g,'/blob_features');
                    % filter green channel by blob size and intensity
                    trajData_g.filtered = (blobFeats_g.area*pixelsize^2<=maxBlobSize_g)&...
                        (blobFeats_g.intensity_mean>=intensityThresholds_g(numCtr));
                end
                if ~strcmp(wormnum,'1W') % if it is multiworm data, we need to filter for worms in clusters
                    num_close_neighbours_rg = h5read(filename,'/num_close_neighbours_rg');
                    inCluster = num_close_neighbours_rg>=minNumNeighbours;
                    features.filtered = ismember(features.skeleton_id+1,find(inCluster));
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
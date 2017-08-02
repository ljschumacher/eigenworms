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

minNumNeighbrs = 3;
minNeighbrDist = 1500;
strains = {'N2','npr1'};
wormnums = {'40','HD'};
analysisTypes = {'loneWorms','inCluster','smallCluster','leaveCluster'};
pixelsize = 100/19.5; % 100 microns are 19.5 pixels

for numCtr = 1:length(wormnums)
    wormnum = wormnums{numCtr};
    figure
    for strainCtr = 1:length(strains)
        strain = strains{strainCtr};
        %% load data
        filenames_r = importdata(['datalists/' strains{strainCtr} '_' wormnum '_r_list.txt']);
        numFiles = length(filenames_r);
        for analysisType = analysisTypes
            peakProjections1 = cell(numFiles,1);
            peakProjections2 = cell(numFiles,1);
            peakProjections3 = cell(numFiles,1);
            peakProjections4 = cell(numFiles,1);
            for fileCtr = 1:numFiles
                filename_r = filenames_r{fileCtr};
                if exist(strrep(filename_r,'skeletons','features'),'file')&&exist(filename_r,'file')
                    features = h5read(strrep(filename_r,'skeletons','features'),'/features_timeseries');
                    trajData = h5read(filename_r,'/trajectories_data');
                    frameRate = h5readatt(filename_r,'/plate_worms','expected_fps');
                    % filter for in-cluster etc
                    min_neighbr_dist = h5read(filename_r,'/min_neighbr_dist');
                    num_close_neighbrs = h5read(filename_r,'/num_close_neighbrs');
                    neighbr_dist = h5read(filename_r,'/neighbr_distances');
                    loneWorms = min_neighbr_dist>=minNeighbrDist;
                    inCluster = num_close_neighbrs>=minNumNeighbrs;
                    smallCluster = (num_close_neighbrs==2 & neighbr_dist(:,3)>=minNeighbrDist)...
                        |(num_close_neighbrs==3 & neighbr_dist(:,4)>=minNeighbrDist)...
                        |(num_close_neighbrs==4 & neighbr_dist(:,5)>=minNeighbrDist);
                    leaveCluster = [false; inCluster(1:end-1)&~inCluster(2:end)];
                    % add 5 more second after each leaving event
                    leaveClusterExtended = unique(find(leaveCluster) + [0:5*double(frameRate)]);
                    leaveClusterExtended = leaveClusterExtended(leaveClusterExtended<numel(leaveCluster)); % exclude frames beyond highest frame number
                    leaveCluster(leaveClusterExtended) = true;
                    leaveCluster(inCluster) = false; % exclude times when worm moved back
                    leaveCluster(loneWorms)=false; % exclude worms that have become lone worm
                    switch analysisType{1}
                        case 'loneWorms'
                            features.filtered = ismember(features.skeleton_id+1,find(loneWorms));
                        case 'inCluster'
                            features.filtered = ismember(features.skeleton_id+1,find(inCluster));
                        case 'smallCluster'
                            features.filtered = ismember(features.skeleton_id+1,find(smallCluster));
                        case 'leaveCluster'
                            features.filtered = ismember(features.skeleton_id+1,find(leaveCluster));
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
                        warning(['All worms filtered out for ' filename_r ])
                    end
                else
                    warning(['Not all necessary tracking results present for ' filename_r ])
                end
            end
            % pool data from multiple recordings
            peakProjections1 = vertcat(peakProjections1{:});
            peakProjections2 = vertcat(peakProjections2{:});
            peakProjections3 = vertcat(peakProjections3{:});
            peakProjections4 = vertcat(peakProjections4{:});
            %% plot data
            subplot(4,2,strainCtr)
            if strcmp(analysisType{1},'loneWorms')
                hold on, title(strain,'FontWeight','normal')
                ylabel('P'), ylim([0 0.28]), yticks(get(gca,'ylim')), xlim([-10 10]), xlabel('a^*_1'), box on
            end
            histogram(peakProjections1,'Normalization','pdf','DisplayStyle','stairs')
            
            subplot(4,2,strainCtr+2)
            if strcmp(analysisType{1},'loneWorms')
                hold on, ylabel('P'), ylim([0 0.28]), yticks(get(gca,'ylim')), xlim([-10 10]), xlabel('a^*_2'), box on
            end
            histogram(peakProjections2,'Normalization','pdf','DisplayStyle','stairs')
            
            subplot(4,2,strainCtr+4)
            if strcmp(analysisType{1},'loneWorms')
                hold on, ylabel('P'), ylim([0 0.28]), yticks(get(gca,'ylim')), xlim([-10 10]), xlabel('a^*_3'), box on
            end
            histogram(peakProjections3,'Normalization','pdf','DisplayStyle','stairs')
            
            subplot(4,2,strainCtr+6)
            if strcmp(analysisType{1},'loneWorms')
                hold on, ylabel('P'), ylim([0 0.28]), yticks(get(gca,'ylim')), xlim([-10 10]), xlabel('a^*_4'), box on
            end
            histogram(peakProjections4,'Normalization','pdf','DisplayStyle','stairs')
        end
    end
    %% format and export figure
    set(gcf,'Position',get(gcf,'Position').*[1 1 1 1.5]) % make figure higher
    legend(analysisTypes)
    figFileName = ['figures/peakAmplitudeEigenwormProjectionDistribution_' wormnum '.eps'];
    exportfig(gcf,figFileName,'Color','rgb')
    system(['epstopdf ' figFileName]);
    system(['rm ' figFileName]);
end
%tilefigs()
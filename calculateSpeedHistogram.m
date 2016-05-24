close all
clear

areas = [];

dBin = 0.05;
bins = 0:dBin:1;
% select data set by strain - N2, HW, NP
for strain = {'N2', 'HW', 'NP'}
    S = strain{:};
    % select data set by number of worms - 1, 5, 15, 25, 40
    for N = [1 5 15 25 40]
        % %         close all
        
        % load file name descriptor - taken from Camille's Recording_LOG.xls
        load(['recordingsLog/strain' S 'worms' num2str(N) '.mat']);
        nFiles = length(filenames);
        % not all results may be present, so check how many
        dsCtr = 1;
        datasets = {};
        for fileCtr=1:nFiles
            % find full path to folder
            file = rdir(['/data1/linus/Recordings/Results/*/' ...
                filenames{fileCtr}(1:end-5) '_skeletons.hdf5']);
            if ~isempty(file)
                datasets{dsCtr} = file.name;
                dsCtr = dsCtr + 1;
            end
        end
        
        if ~isempty(datasets)
            nDatasets = length(datasets);
            
            for dsCtr=1:nDatasets
                filename = datasets{dsCtr};
                % load trajectory data
                trajectoryData = h5read(filename,'/trajectories_data');
                hasSkel = trajectoryData.has_skeleton==1;
                framesPerWorm = histcounts(trajectoryData.worm_index_joined,max(trajectoryData.worm_index_joined));
                frequentWorms = find(framesPerWorm>=250);
                frequentFilter = ismember(trajectoryData.worm_index_joined,frequentWorms);
                % select frames with a certain area
                dABin = 25;
                Abins = 0:dABin:1500;
                counts = histcounts(trajectoryData.area(hasSkel&frequentFilter),Abins,'normalization','pdf');
                [~, locs, widths, proms] = findpeaks(counts,Abins(2:end)-dABin/2,...
                    'MinPeakDistance',50,'MinPeakWidth',50);
                [~, mostProm] = max(proms); % find most prominent peak
                areaFilter = trajectoryData.area>=(locs(mostProm) - widths(mostProm))&...
                    trajectoryData.area<=(locs(mostProm) + widths(mostProm));
                % find individual worm trajectories
                wormIDs = unique(trajectoryData.worm_index_joined(hasSkel&frequentFilter&areaFilter))';
                figure, hold on
                % go through worms and calculate speeds
                for wormCtr=wormIDs
                    wormIdcs = trajectoryData.worm_index_joined==wormCtr;
                    wormDx = diff(trajectoryData.coord_x(wormIdcs));
                    wormDy = diff(trajectoryData.coord_y(wormIdcs));
                    wormDf = diff(trajectoryData.frame_number(wormIdcs));
                    wormDisplacement = sqrt(wormDx.^2 + wormDy.^2)./single(wormDf);
                    histogram(wormDisplacement,bins)
                    % if isfield(trajData,'worm_label')
                    %                 wormLabels = unique(trajData.worm_label(hasSkel&frequentFilter&areaFilter))';
                    %                 figure, hold on
                    %                 % go through worms and calculate speeds
                    %                 for wormLbl=wormLabels
                    %                     wormIdcs = trajData.worm_label==wormLbl;
                    %                     wormDx = diff(trajData.coord_x(wormIdcs));
                    %                     wormDy = diff(trajData.coord_y(wormIdcs));
                    %                     wormDisplacement = sqrt(wormDx.^2 + wormDy.^2);
                    %                     histogram(wormDisplacement,bins)
                    %                 end
                    % end
                    xlim([0 1])
                    title([S ' N=' num2str(N)])
                    set(gcf,'name',filename)
                end
                pause(0.001)
            end
        else
            display(['No datasets found for strain=' S ', worms=' num2str(N) ])
        end
    end
end
% histogram(areas)
tilefigs%([5 6])
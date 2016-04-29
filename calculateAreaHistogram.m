close all
clear

areas = [];

dBin = 25;
bins = 0:dBin:1500;
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
                framesFilter = ismember(trajectoryData.worm_index_joined,frequentWorms);
                figure
                h = histogram(trajectoryData.area(hasSkel&framesFilter),bins,'normalization','pdf');
                title([S ' N=' num2str(N)])
                set(gcf,'name',filename)
                xlim([0 2500])
                hold on
                [peaks, locs, widths, proms] = findpeaks(h.Values,double(h.BinEdges(2:end))-dBin/2,...
                    'MinPeakDistance',50,'MinPeakWidth',50);
                plot(locs,peaks,'v','LineWidth',2)
                % find most prominent peak
                [~, mostProm] = max(proms);
                plot(locs(mostProm),peaks(mostProm),'ro','LineWidth',2,'MarkerSize',20)
                stairs([locs - widths; locs - widths; locs + widths],...
                    [0; 1; 0]*peaks,'LineWidth',2)
            end
        else
            display(['No datasets found for strain=' S ', worms=' num2str(N) ])
        end
    end
end
% histogram(areas)
tilefigs([5 6])
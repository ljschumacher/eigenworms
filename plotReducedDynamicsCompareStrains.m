% combine projected amplitudes in space of first two eigenworms

% issues/to-do:
% - should we normalise the variance of each projection to 1, as in the
% original eigenworm paper?
close all
clear

% figure export options
exportOptions = struct('Color','rgb','Renderer','painters');

strains = {'N2', 'HW', 'NP'};
nStrains = length(strains);
for N = [1 5 15 25 40]
    % loop through different strains
    for strainCtr = 1:nStrains
        S = strains{strainCtr};
        % load eigenworm data
        file = rdir(['results/' S '_' num2str(N) 'worms_eigenData.mat']);
        if ~isempty(file)
            eigProjectionFig = figure;
            % load eigenworm analysis result
            load(file.name)
            % normalise to unit variance
            eigenProjections = zscore(eigenProjections);
            % plot projected amplitudes
            subplot(1,3,1)
            histogram2(eigenProjections(:,1),eigenProjections(:,2),...
                'Normalization','pdf','DisplayStyle','tile','EdgeColor','none','FaceColor','flat')
            title(['strain ' S])
            xlabel('a_1'), ylabel('a_2')
            axis image
            grid off
            xlim([-2 2]), ylim([-2 2])
            % PC1-3
            subplot(1,3,2)
            histogram2(eigenProjections(:,1),eigenProjections(:,3),...
                'Normalization','pdf','DisplayStyle','tile','EdgeColor','none','FaceColor','flat')
            xlabel('a_1'), ylabel('a_3')
            axis image
            grid off
            xlim([-2 2]), ylim([-2 2])
            title(['N = ' num2str(N) ' worms'])
            % PC2-3
            subplot(1,3,3)
            histogram2(eigenProjections(:,2),eigenProjections(:,3),...
                'Normalization','pdf','DisplayStyle','tile','EdgeColor','none','FaceColor','flat')
            xlabel('a_2'), ylabel('a_3')
            axis image
            grid off
            xlim([-2 2]), ylim([-2 2])
            title([num2str(size(eigenProjections,1)/25/3600,2) ' worm-hours'])
            % annotate and save figure
            set(eigProjectionFig, 'name', ['projected amplitudes for N=' num2str(N) ' worms'])
            figName = ['figures/' S '_' num2str(N) 'worms_projections.eps'];
            exportfig(eigProjectionFig,figName,exportOptions)
            system(['epstopdf ' figName]);
            system(['rm ' figName]);
%             close(eigProjectionFig)
            clear eigenProjections
        else
            display(['No data for strain=' S ', worms=' num2str(N)])
        end
    end
end
tilefigs([4 3])
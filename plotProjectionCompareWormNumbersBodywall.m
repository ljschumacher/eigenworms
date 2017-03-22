% combine projected amplitudes in space of first two eigenworms

% issues/to-do:
% -
close all
clear

% figure export options
exportOptions = struct('Color','rgb');

nEigenworms = 4;

strains = {'N2', 'npr1'};
nStrains = length(strains);
wormnums = {'HD','40','1W'};
for strainCtr = 1:nStrains
    S = strains{strainCtr};
    eigProjectionFig = figure;
    plotHandles = NaN(size(wormnums));
    for numCtr = 1:length(wormnums)
        N = wormnums{numCtr};
        % load eigenworm data
        file = rdir(['results/eigenData_' S '_' N '_bodywall.mat']);
        if ~isempty(file)
            % load eigenworm analysis result
            load(file.name,'masterProjections');
%             % change the order of first and third reference eigenworm
%             masterProjections = masterProjections(:,[3 2 1 4:end]);
            % normalise to unit variance
            % %             masterProjections = zscore(masterProjections);
            % plot projected amplitudes
            for cmpCtr = 1:nEigenworms
                subplot(ceil(nEigenworms/2),2,cmpCtr)
                plotHandles(numCtr) =...
                    histogram(masterProjections(:,cmpCtr),...
                    'Normalization','pdf','DisplayStyle','stairs');
                if numCtr == 1
                    hold on
                    ax = gca;
                    ax.YLabel.String = 'P';
                    ax.XLabel.String =  ['a_' num2str(cmpCtr)];
                    ax.XLim = [-10 10];
                end
            end
            clear masterProjections
        else
            display(['No data for strain=' S ', worms=' num2str(N)])
        end
    end
    % annotate and save figure
    legend(plotHandles,wormnums)
    set(eigProjectionFig, 'name', ['projected amplitudes ' S])
    figName = ['projections_' S '_CompareWormNumbers'];
    exportfig(eigProjectionFig,['figures/' figName '.eps'],exportOptions)
    system(['epstopdf figures/' figName '.eps']);
    system(['rm figures/' figName '.eps']);
    %             close(eigProjectionFig)
end
tilefigs([2 3])
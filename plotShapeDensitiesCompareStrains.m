% combine projected amplitudes in space of first two eigenworms
% issues / to-do:
% - should the contour levels be at some cumulative probability, i.e.
% whatever captures 90% or so of poses?
close all
clear

% figure export options
exportOptions = struct('Color','rgb');

nComponents = 4;

strains = {'N2', 'HW', 'NP'};
nStrains = length(strains);
plotColors = lines(nStrains);
for N = [1 5 15 25 40]
    eigProjectionFig = figure;
    hold on
    % loop through different strains
    for strainCtr = 1:nStrains
        S = strains{strainCtr};
        % load eigenworm data
        file = rdir(['results/' S '_' num2str(N) 'worms_eigenData.mat']);
        if ~isempty(file)
            % load eigenworm analysis result
            load(file.name,'masterProjections')
            % normalise to unit variance
            masterProjections = zscore(masterProjections);
            % calculate histogram
            h = histogram2(masterProjections(:,2),masterProjections(:,3),...
                'Normalization','pdf','Visible','off');
            [~, lineHandles(strainCtr)] = contour(h.XBinEdges(2:end)-h.BinWidth(1)/2,...
                h.YBinEdges(2:end)-h.BinWidth(2)/2,h.Values',...
                max(max(h.Values/2)).*[1 1],'Color',plotColors(strainCtr,:),...
                'LineWidth',2);
            clear masterProjections
        else
            display(['No data for strain=' S ', worms=' num2str(N)])
        end
    end
    % annotate and save figure
    lh = legend(lineHandles,strains);
    lh.Title.String = 'strain';
    title(['N = ' num2str(N) ' worms'])
    axis equal
    ylim([-2 2])
    xlim([-2 2])
    box on
    xlabel('a_1'), ylabel('a_2')
    set(eigProjectionFig, 'name', ['projected amplitudes for N=' num2str(N) ' worms'])
    figName = ['figures/projections_' num2str(N) 'worms.eps'];
    exportfig(eigProjectionFig,figName,exportOptions)
    system(['epstopdf ' figName]);
    system(['rm ' figName]);
    %             close(eigProjectionFig)
end
tilefigs()
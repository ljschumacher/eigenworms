% combine projected amplitudes in space of first two eigenworms

% issues/to-do:
% -
close all
clear

% figure export options
exportOptions = struct('Color','rgb');

strains = {'npr1','N2'};
nStrains = length(strains);
wormnums = {'40'}%{'HD','40'};
plotColors = lines(nStrains);

for numCtr = 1:length(wormnums)
    N = wormnums{numCtr};
    varFig = figure;
    for strainCtr = 1:nStrains
        S = strains{strainCtr};
        % load eigenworm data
        file = rdir(['results/eigenData_' S '_' N '_bodywall_loneWorms_joining.mat']);
        % load eigenworm analysis result
        load(file.name,'masterVarExplained');
        plot(0:6,cumsum([0; masterVarExplained]),'--o','Color',plotColors(strainCtr,:))
        hold on
        % load eigenworm data
        file = rdir(['results/eigenData_' S '_' N '_bodywall_inCluster_joining.mat']);
        % load eigenworm analysis result
        load(file.name,'masterVarExplained');
        plot(0:6,cumsum([0; masterVarExplained]),'--x','Color',plotColors(strainCtr,:))
    end
    ax = gca;
    ax.YLabel.String = '% variance explained';
    ax.XLabel.String = 'eigenworm';
    % annotate and save figure
    legend({'npr-1 lone','npr-1 cluster','N2 lone','N2 cluster'},'Location','SouthEast')
    title(['variance explained for lone (-o) and in cluster (-x) eigenworms'],'FontWeight','normal')
    figName = ['figures/varExplained_CompareClusterStatus.eps'];
    exportfig(varFig,figName,exportOptions)
    system(['epstopdf ' figName]);
    system(['rm ' figName]);
end

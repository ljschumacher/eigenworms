% combine eigenworm plots from saved eigenworms for different strains

% issues/to-do
% - to plot variance explained, we actually need to save all eigenvalues,
% or at least their sum. all eigenvalues might also be useful for
% statistical tests on their distribution.
close all
clear

nEigenworms = 6;
N = 40;
strains = {'N2', 'HW', 'NP'};

eigWormFig = figure;
for eigCtr = 1:nEigenworms
    subplot(ceil(nEigenworms/2),2,eigCtr)
    hold on
    box on
    title(num2str(eigCtr))
end

% loop through different strains
for S = strains
    % load eigenworm data
    file = rdir(['results/' S{1} '_' num2str(N) 'worms_eigenData.mat']);
    if ~isempty(file)
        % load eigenworm analysis result
        load(file.name)
        % plot first few eigenworms
        for eigCtr = 1:nEigenworms
            plot(eigWormFig.Children(nEigenworms - eigCtr + 1),...
                eigenWorms(eigCtr,:),'LineWidth',2)
        end
    else
        display(['No data for strain=' S{1} ', worms=' num2str(N)])
    end
end
% annotate and save figure
legH = legend(eigWormFig.Children(1),strains);
legH.Title.String = 'strain';
set(eigWormFig, 'name', ['eigenworms for N=' num2str(N) ' worms'])
savefig(eigWormFig,['figures/' num2str(N) 'worms_eigenworms_compareStrains.fig'],'compact')

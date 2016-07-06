% combine eigenworm plots from saved eigenworms for different worm numbers

% issues/to-do
% - to plot variance explained, we actually need to save all eigenvalues,
% or at least their sum. all eigenvalues might also be useful for
% statistical tests on their distribution.
close all
clear

% figure export options
exportOptions = struct('Color','rgb');

nEigenworms = 6;
wormNumbers = [1 5 15 25 40];
nWormNumbers = length(wormNumbers);
plotColors= lines(nWormNumbers);

% for reference, eigenworms from Brown et al. 2013
reference = load('../Motif_Analysis/eigenWorms.mat');
% change the order of first and third masterworm to match ours
reference.eigenWorms = reference.eigenWorms([3 2 1 4:end], :);

% loop through different strains
for strain = {'N2', 'HW', 'NP'}
    S = strain{:};
    eigWormFig = figure;
    for eigCtr = 1:nEigenworms
        subplot(ceil(nEigenworms/2),2,eigCtr)
        hold on
        box on
        title(num2str(eigCtr))
    end
    % loop through worm numbers - 1, 5, 15, 25, 40
    missingN = false(1,nWormNumbers);
    for numCtr = 1:nWormNumbers;
        N = wormNumbers(numCtr);
        % load eigenworm data
        file = rdir(['results/' S '_' num2str(N) 'worms_eigenData.mat']);
        if ~isempty(file)
            % load eigenworm analysis result
            load(file.name,'eigenWorms')
            % plot first few eigenworms
            for eigCtr = 1:nEigenworms
                % check if we need to flip eigenworm upside-down
                if eigenWorms(eigCtr,:)*reference.eigenWorms(eigCtr,:)'<0
                    eigenWorms(eigCtr,:) = -eigenWorms(eigCtr,:);
                end
                plot(eigWormFig.Children(nEigenworms - eigCtr + 1),...
                    eigenWorms(eigCtr,:),'Color', plotColors(numCtr,:),...
                    'LineWidth',2)
            end
        else
            display(['No data for strain=' S ', worms=' num2str(N)])
            missingN(numCtr) = true;
        end
    end
    % plot reference eigenworms from Brown et al. 2013
    for eigCtr = 1:nEigenworms
        subplot(ceil(nEigenworms/2),2,eigCtr)
        plot(reference.eigenWorms(eigCtr,:),'k--')
    end
    % annotate and save figure
    legH = legend(eigWormFig.Children(1),num2str(wormNumbers(~missingN)'));
    %     legH.Title.String = 'N worms';
    set(eigWormFig, 'name', ['Strain ' S ' eigenworms'])
    figFileName = ['figures/' S '_eigenworms_compareN.eps'];
    exportfig(eigWormFig,figFileName,exportOptions)
    system(['epstopdf ' figFileName]);
    system(['rm ' figFileName]);
end
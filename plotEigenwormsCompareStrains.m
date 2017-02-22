% combine eigenworm plots from saved eigenworms for different strains

% issues/to-do

close all
clear

% figure export options
exportOptions = struct('Color','rgb');

% for reference, eigenworms from Brown et al. 2013
reference = load('../Motif_Analysis/eigenWorms.mat');
% change the order of first and third masterworm to match ours
reference.eigenWorms = reference.eigenWorms([3 2 1 4:end], :);

nEigenworms = 6;
for N = [1 5 15 25 40]
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
            load(file.name,'eigenWorms')
            % plot first few eigenworms
            for eigCtr = 1:nEigenworms
                % check if we need to flip eigenworm upside-down
                if eigenWorms(eigCtr,:)*reference.eigenWorms(eigCtr,:)'<0
                    eigenWorms(eigCtr,:) = -eigenWorms(eigCtr,:);
                end
                plot(eigWormFig.Children(nEigenworms - eigCtr + 1),...
                    eigenWorms(eigCtr,:),'LineWidth',2)
            end
        else
            display(['No data for strain=' S{1} ', worms=' num2str(N)])
        end
    end
    % plot reference eigenworms from Brown et al. 2013
    for eigCtr = 1:nEigenworms
        subplot(ceil(nEigenworms/2),2,eigCtr)
        plot(reference.eigenWorms(eigCtr,:),'k--')
    end
    % annotate and save figure
    legH = legend(eigWormFig.Children(1),strains);
    %     legH.Title.String = 'strain';
    set(eigWormFig, 'name', ['eigenworms for N=' num2str(N) ' worms'])
    figFileName = ['figures/' num2str(N) 'worms_eigenworms_compareStrains.eps'];
    exportfig(eigWormFig,figFileName,exportOptions)
    system(['epstopdf ' figFileName]);
    system(['rm ' figFileName]);
end
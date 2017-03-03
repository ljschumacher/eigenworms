% combine eigenworm plots from saved eigenworms for different worm numbers
% plotting eigenworms in real space rather than in angles
% plotting variation of eigenworms

close all
clear

% figure export options
exportOptions = struct('Color','rgb');

nEigenworms = 4;
strains = {'N2', 'npr1'};
nStrains = length(strains);
wormnums = {'HD','40','1W'};
nVariations = 5;
weights = linspace(-1,1,nVariations)';

% for reference, eigenworms from Brown et al. 2013
reference = load('../Motif_Analysis/eigenWorms.mat');
% change the order of first and third masterworm to match ours
% reference.eigenWorms = reference.eigenWorms([3 2 1 4:end], :);

%% change the amount being varied to 2 sigma of the projection (own)
                %sigma2 = 2*std(masterProjections(:,cmpCtr));
%% also use this to plot variation on projections of master eigenworms                
for strainCtr = 1:nStrains
    S = strains{strainCtr};
    eigWormFig = figure;
    for numCtr = 1:length(wormnums)
        N = wormnums{numCtr};
        % load eigenworm data
        file = rdir(['results/eigenData_' S '_' N '_bodywall.mat']);
        if ~isempty(file)
            % load eigenworm analysis result
            load(file.name,'eigenWorms')
            % plot first few eigenworms in real space
            for eigCtr = 1:nEigenworms
                subplot(nEigenworms,length(wormnums),(eigCtr-1)*length(wormnums)+numCtr), hold on
                set(gca,'ColorOrder',parula(nVariations),'Color','none')
                [x, y] = angles2xy(weights*eigenWorms(eigCtr,:));
                % check if we need to flip eigenworm upside-down
                if eigenWorms(eigCtr,:)*reference.eigenWorms(eigCtr,:)'<0
                    eigenWorms(eigCtr,:) = -eigenWorms(eigCtr,:);
                end
                plot(x',y','LineWidth',2)
                if eigCtr==1
                    title(N,'FontWeight','normal')
                end
                if numCtr==1
                    ylabel(num2str(eigCtr))
                end
                xlim([0 48])
                axis equal
                set(gca,'XTick',[],'YTick',[])
            end
        else
            display(['No data for strain=' S ', worms=' num2str(N)])
        end
    end
    % annotate and save figure
    set(eigWormFig, 'name', ['Strain ' S ' ' N ' eigenworms'])
    figFileName = ['figures/' S '_' N '_eigenworms_variation.eps'];
    exportfig(eigWormFig,figFileName,exportOptions)
    system(['epstopdf ' figFileName]);
    system(['rm ' figFileName]);
end
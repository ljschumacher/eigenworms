% plot the correlations of tangent vectors (or average cos theta) also a
% worm backbone

close all
clear

% set parameters for correlation analysis
maxLag = 40;

strains = {'npr-1'};

for strain = strains
    % find all files in the directory tree
    files = rdir(['singleWorm/' strain{:} '/**/*.mat']);
    numWorms = size(files,1);
    
    tanCorr = NaN(numWorms,maxLag + 1);
    
    for wormCtr = 1:numWorms
        % load single worm data
        load(files(wormCtr).name)
        
        % clean worm posture data from NaNs
        cleanX = worm.posture.skeleton.x(:,~any(isnan(worm.posture.skeleton.x),1));
        cleanY = worm.posture.skeleton.y(:,~any(isnan(worm.posture.skeleton.x),1));
        numSamples = size(cleanX,2);
        
        swCorr = NaN(numSamples,maxLag + 1);
        
        % convert positions into components of tangent vectors
        dX = diff(cleanX,1,1);
        dY = diff(cleanY,1,1);
        
        % loop through each sample and calculate the dot product of angles - can
        % probably speed this up by calculating the dot-product elementwise for all
        % samples at the same time (for each lag)
        for sampleCtr = 1:numSamples
            swCorr(sampleCtr,:) = vectorAutoCorrelation([dX(:,sampleCtr), dY(:,sampleCtr)]',maxLag,0,1,1);
        end
        % take the average over all samples from this worm
        tanCorr(wormCtr,:) = mean(swCorr,1);
    end
    
    % plot results
    plot(0:maxLag,tanCorr)
    ylim([0 1])
    
    % annotate and save figure
    ylabel('$\langle\cos\theta\rangle$','Interpreter','LaTeX')
    xlabel('$\Delta s$','Interpreter','LaTeX')
    xlim([0 maxLag])
    figName = [strain{:} '_wormPersistence.eps'];
    exportfig(gcf,figName,'color','rgb')
    system(['epstopdf ' figName]);
    system(['rm ' figName]);
end

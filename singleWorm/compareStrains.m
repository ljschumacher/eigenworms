% plot mean eigen-projections in multiple dimensions to compare strains
close all
clear

files = dir('*.mat');

nStrains = size(files,1);

meanProj = NaN(nStrains,6);
errProj = NaN(size(meanProj));

legendStrings = cell(nStrains,1);

for strainCtr = 1:nStrains
    load(files(strainCtr).name)
    legendStrings{strainCtr} = strrep(strrep(strrep(files(strainCtr).name,...
        '_',' '),'histogram',''),'.mat','');
    for eigCtr = 1:6
        meanProj(strainCtr,eigCtr) = ...
            worm.posture.eigenProjection(eigCtr).histogram.sets.mean.abs;
        errProj(strainCtr,eigCtr) = ...
            worm.posture.eigenProjection(eigCtr).histogram.sets.stdDev.abs...
            /sqrt(worm.posture.eigenProjection(eigCtr).histogram.sets.samples - 1);
    end
end

errorbar(meanProj',errProj')
legend(legendStrings)
ylabel('\langle|eigen projection|\rangle')
xlabel('eigenworm number')
xlim([1 6])
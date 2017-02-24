function [ ] = plotSkelLengthDist(skelData,pixelsize,minSkelLength,filename)


skelLengthFig = figure;
histogram(sum(sqrt(sum((diff(skelData,1,2)*pixelsize).^2))),...
    'Normalization','Probability','EdgeColor','none')
hold on
plot([minSkelLength minSkelLength],[0 0.1],'r--')
xlabel('skeleton length (\mu m)'), ylabel('P')
figName = strrep(strrep(filename(end-31:end),'/',''),'_',' ');
title(figName,'Fontweight','normal')

% figure export
figFileName = ['figures/diagnostics/' strrep(strrep(figName,' ','_'),'s.hdf5','') 'lengths.eps'];
exportfig(skelLengthFig,figFileName,'Color','rgb')
system(['epstopdf ' figFileName]);
system(['rm ' figFileName]);

close(skelLengthFig)

end


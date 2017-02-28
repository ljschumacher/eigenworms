function [ ] = plotSkelLengthDist(skelData,pixelsize,minSkelLength,filename)


skelLengthFig = figure;
skelLengths = sum(sqrt(sum((diff(skelData,1,2)*pixelsize).^2)));
histogram(skelLengths,'Normalization','Probability','EdgeColor','none')
hold on
plot([minSkelLength minSkelLength],[0 0.1],'r--')
xlabel('skeleton length (\mu m)'), ylabel('P')
figName = strrep(filename,'_',' ');
title(figName,'Fontweight','normal')

% plot filtered skeletons on top
axpos = get(gca,'Position');
newpos = axpos;
newpos(3) = 0.33*axpos(3); % third width
overlayRejected = axes('Position',newpos,'Visible','off');
plotSkeletons(skelData(:,:,skelLengths<minSkelLength),9,'r');
newpos(1) = axpos(1) + 0.67*axpos(3); % start at 2/3 width
overlayAccepted = axes('Position',newpos,'Visible','off');
plotSkeletons(skelData(:,:,skelLengths>minSkelLength),9,'k');
% figure export
figFileName = ['figures/diagnostics/' strrep(strrep(figName,' ','_'),'s.hdf5','') 'lengths.eps'];
exportfig(skelLengthFig,figFileName,'Color','rgb')
system(['epstopdf ' figFileName]);
system(['rm ' figFileName]);

close(skelLengthFig)

end


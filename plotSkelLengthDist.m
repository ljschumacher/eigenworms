function [ ] = plotSkelLengthDist(skelData,pixelsize,minSkelLength,maxSkelLength,filename)


skelLengthFig = figure;
skelLengths = sum(sqrt(sum((diff(skelData,1,2)*pixelsize).^2)));
histogram(skelLengths,'Normalization','pdf','EdgeColor','none')
hold on
plot(minSkelLength,[0 0.02],'r--')
plot(maxSkelLength,[0 0.02],'r--')
xlabel('skeleton length (\mu m)'), ylabel('P')
figName = strrep(filename,'_',' ');
title(figName,'Fontweight','normal')

% plot filtered skeletons on top
axpos = get(gca,'Position');
rejected = skelLengths<minSkelLength;
newpos = axpos;
newpos(3) = 0.33*axpos(3); % third width
if any(rejected)
    overlayRejected = axes('Position',newpos,'Visible','off');
    plotSkeletons(skelData(:,:,rejected),min(9,nnz(rejected)),'r');
end
accepted = skelLengths>minSkelLength;
if any(accepted)
    newpos(1) = axpos(1) + 0.67*axpos(3); % start at 2/3 width
    overlayAccepted = axes('Position',newpos,'Visible','off');
    plotSkeletons(skelData(:,:,skelLengths>minSkelLength),min(9,nnz(accepted)),'k');
else
    1;
end
% figure export
figFileName = ['figures/diagnostics/' strrep(strrep(figName,' ','_'),'s.hdf5','') 'lengths.eps'];
exportfig(skelLengthFig,figFileName,'Color','rgb')
system(['epstopdf ' figFileName]);
system(['rm ' figFileName]);

close(skelLengthFig)

end


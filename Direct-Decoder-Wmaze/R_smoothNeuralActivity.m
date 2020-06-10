function [ySum,yMean, yCells ]= R_smoothNeuralActivity(Spike_Binary, TimeWindow)
[nt,nc]=size(Spike_Binary);
TimeWindow=repmat(TimeWindow,nc,1);
yCells=conv2(Spike_Binary,TimeWindow','same');
yMean=mean(yCells,2);
ySum=sum(yCells,2);
%X=conv(X,ones(LenWind,1)/LenWind,'same');
%X=downsample(X,LenWind);
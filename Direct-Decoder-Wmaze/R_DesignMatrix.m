function Xdsgn = R_DesignMatrix(X,NumFilters)
%% Create design matrix
% Build design matrix for the stimulus
[times,channels]=size(X);
PadX = [zeros(NumFilters-1,channels); X];
for i=1:channels
Xdsgn(:,(i-1)*NumFilters+1:i*NumFilters) = hankel(PadX(1:end-NumFilters+1,i),...
     PadX(end-NumFilters+1:end,i));

end

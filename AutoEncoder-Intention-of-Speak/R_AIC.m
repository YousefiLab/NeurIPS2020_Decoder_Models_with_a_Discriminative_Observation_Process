function [AIC,BIC]=R_AIC(ratepred_pGLM,b, yMean)
Ysp=zeros(size(yMean));
Ysp(find(yMean ~= 0))=1;
%% AIC measure for goodness of fit measure[Stone, 1979; Akaike, 1974]
LL_histGLM = Ysp'*log(ratepred_pGLM) - sum(ratepred_pGLM);
AIC = -2*LL_histGLM + 2*(length(b));
BIC = -2*LL_histGLM + log(sum(Ysp))*(length(b));

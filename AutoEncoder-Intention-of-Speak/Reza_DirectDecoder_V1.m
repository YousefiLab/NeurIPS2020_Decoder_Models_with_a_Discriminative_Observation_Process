% close all
clear all
% clc
%% rebuild neural feature matrix using original data window (train+train)
%Reza S is neureal feature, X is cognetive state, Z is behaviour
load('Data3.mat');
load('Apr12.mat');

%% rebuild neural feature matrix using original data window (train+train)
%Reza S is neureal feature, X is cognetive state, Z is behaviour

nTimeValues = numel(Data.train.t);
nStateValues = numel(xvals);
P = nan(nTimeValues,nStateValues);
Pzx = nan(nTimeValues,nStateValues);
PSmooth = nan(nTimeValues,nStateValues);
Pxxkn1YT=nan(nStateValues,nStateValues);
xkHat = nan(nTimeValues,1);
xs = nan(nTimeValues,1);

% Pxx = nan(nStateValues);
Pxxc = nan(nStateValues);



%%% Pxx
for ii = 1:nStateValues % total number of states 
    mxx = Par.XX.Ak*xvals(ii) + Par.XX.Bk;
%     Pxx(ii,:) = normpdf(xvals,mxx,Par.XX.sigma);
    Pxxc(ii,:) = normpdf(xvals,mxx,Par.XX.sigma);
end
for ii = 1:nStateValues
%     Pxx(:,ii) = Pxx(:,ii)/sum(Pxx(:,ii));
    Pxxc(:,ii) = Pxxc(:,ii)/sum(Pxxc(:,ii));
end
%%% Pxs
Pxs_prv = ones(size(Pxxc(1,:)))/length(Pxxc(1,:));

for ii = 1:nTimeValues  
      
    if ii == 1
      
      Px = Pxxc * Pxs_prv'; % integral(p(xk|xk-1) * Pxs)
      Px = Px'/sum(Px');
      P(ii,:) = Px;
       
    else
       
%          ii
         Px = Pxxc * Pxs_prv'; % integral(p(xk|xk-1) * Pxs)
         Px = Px'/sum(Px');
        %%%%%%%% predict step %%%%%%%%%%%
        % compute p(xk|s1,...,sk-1) = integral(p(xk|xk-1) * Pprev)
        Pint = Pxxc * P(ii-1,:)';
        
        % compute p(xk|sk,hk) 
       
        xmu = glmval( Par.XSH.w,  Data.train.SH(ii,:), 'identity');
        Pxs = normpdf(xvals,xmu,Par.XSH.sigma);
       
        Pxs_prv = Pxs;
        % compute P = p(xk|sk) * p(xk|s1,...,sk-1) * p(zk | xk) / p(xk)
        
        p = (Pxs .* Pint'  ) ./ Px;
        P(ii,:) =( p/sum(p));
      
       
    end
        
        xkHat(ii) = xvals*P(ii,:)';
        xsHat(ii) = xvals.*xvals*P(ii,:)' - xkHat(ii)^2;
end
%% Smoother
for ii = nTimeValues:-1:1 
    if ii == nTimeValues
        PSmooth(ii,:)=P(ii,:);
        epj(ii)=xkHat(ii);
    else
        Pxkp1Sk=Pxxc * P(ii,:)';
        PSmooth(ii,:)=P(ii,:) .* (Pxxc * (PSmooth(ii+1,:)' ./Pxkp1Sk))';
       
    end
    xkHatS(ii) = xvals*PSmooth(ii,:)';
    xsHatS(ii) = xvals.*xvals*PSmooth(ii,:)' - xkHatS(ii)^2;
end


%%
tt = Data.train.ti-Data.train.ti(1);
tx = Data.train.t-Data.train.t(1);
figure
set(gcf,'color','w')

hold on
% h(1) = plot(tx,Data.train.x,'color',[1 0 1],'linewidth',2);
h(1) = plot(tt,xkHatSAutoEncoder,'color',0*[1 1 1],'linewidth',2);
% he = shadedErrorBar(tt,xkHat,sqrt(xsHat),{'-','color',[0 0 0]},1);
h(2) = plot(tx,xkHatS,'color',[0 1 0],'linewidth',2);
he = shadedErrorBar(tx,xkHatS,sqrt(xsHatS),{'-','color',[0 1 0]},1);
yvalsC = 0*stateBound*ones(size(Data.train.cwords));
h(3) = plot(Data.train.cwords,yvalsC,'or');

yvalsP = 1*stateBound*ones(size(Data.train.pwords));
h(4) = plot(Data.train.pwords,yvalsP,'ob');hold off
hold off
set(gca,'ydir','normal')
ylabel('state','fontsize',22);
xlabel('seconds','fontsize',22);
xlim([min(tt),max(tt)]);
ylim([0*stateBound,1*stateBound]);
title('Estimated State with D-D model','fontsize',22);
legend(h,{'Auto-encoder result','D-D result','Companion','Participant' },'fontsize',14);
% colorbar




%%


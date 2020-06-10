close all
clear all
clc
%% rebuild neural feature matrix using original data window (train+train)
%Reza S is neureal feature, X is cognetive state, Z is behaviour
load('Data3.mat');
dt=.05;
dvd=1/dt*2;
Data.train.ti=0:dt:max(Data.train.t); % 100 ms resolution
% Data.test.ti=0:dt:max(Data.test.t);
hs=1:10;
for bicIndx=1:length(hs)
Data.Hs=hs(bicIndx);
Data.ZRes=0.05;
Data.ZHis=0.1;
% Data.train.x=-1*Data.train.x;
gausWindowLength=30;
gausWindowSigma=4;

gausWindow = fspecial('gaussian', gausWindowLength, gausWindowSigma);
gausWindow=mean(gausWindow);
gausWindow=gausWindow/max(gausWindow);
[ySum,ymean, yCells ]= R_smoothNeuralActivity_v2(real(log(Data.train.S(:,Data.bestIdx))), gausWindow);

Data.train.SH=R_DesignMatrix(yCells,Data.Hs);


%%
stateBound=1;
%% Parameter setting (initialization)
MaxItr=3;
for itr=1:MaxItr
    fprintf('itr %d/%d\n',itr,MaxItr);
    if itr == 1
        Par=parInit(Data,stateBound);
    end

xvals= -stateBound:0.02:stateBound;
nTimeValues = numel(Data.train.ti);
nStateValues = numel(xvals);
P = nan(nTimeValues,nStateValues);
Pzx = nan(nTimeValues,nStateValues);
PSmooth = nan(nTimeValues,nStateValues);
Pxxkn1YT=nan(nStateValues,nStateValues);
xkHat = nan(nTimeValues,1);
xs = nan(nTimeValues,1);
pw = zeros(nTimeValues,1);
qw = zeros(nTimeValues,1);
zk = zeros(nTimeValues,1);
Pxx = nan(nStateValues);
Pxxc = nan(nStateValues);



%%% Pxx
for ii = 1:nStateValues % total number of states 
    mxx = Par.XX.Ak*xvals(ii) + Par.XX.Bk;
    Pxx(ii,:) = normpdf(xvals,mxx,Par.XX.sigma);
    Pxxc(ii,:) = normpdf(xvals,mxx,(sqrt(dvd))*Par.XX.sigma);
end
for ii = 1:nStateValues
    Pxx(:,ii) = Pxx(:,ii)/sum(Pxx(:,ii));
    Pxxc(:,ii) = Pxxc(:,ii)/sum(Pxxc(:,ii));
end
%%% Pxs
Pxs_prv = ones(size(Pxx(1,:)))/length(Pxx(1,:));

for ii = 1:nTimeValues  
      
    if ii == 1
      
      Px = Pxxc * Pxs_prv'; % integral(p(xk|xk-1) * Pxs)
      Px = Px'/sum(Px');
      P(ii,:) = Px;
       
    else
        if (ii-1)/(dvd)-floor((ii-1)/(dvd)) == 0
%          ii
         Px = Pxxc * Pxs_prv'; % integral(p(xk|xk-1) * Pxs)
         Px = Px'/sum(Px');
        %%%%%%%% predict step %%%%%%%%%%%
        % compute p(xk|s1,...,sk-1) = integral(p(xk|xk-1) * Pprev)
        Pint = Pxxc * P(ii-1,:)';
        % compute p(zk|xk)
        [Pzx(ii,:) ,pw(ii),qw(ii),z(ii)]=R_Pzx(Data, Par,xvals,ii);
        % compute p(xk|sk,hk) 
       
        xmu = glmval( Par.XSH.w,  Data.train.SH((ii-1)/(dvd),:), 'identity');
        Pxs = normpdf(xvals,xmu,(sqrt(dvd))*Par.XSH.sigma);
        Pxs=Pxs./sum(Pxs);
        Pxs_prv = Pxs;
        % compute P = p(xk|sk) * p(xk|s1,...,sk-1) * p(zk | xk) / p(xk)
        if sum(Pzx(ii,:)) ~=0
        p = (Pxs .* Pint'.* Pzx(ii,:)  ) ./ Px;
        else
       p = (Pxs .* Pint'  ) ./ Px;
        end
        else
        
        
        %%%%%%%% predict step %%%%%%%%%%%
        % compute p(xk|s1,...,sk-1) = integral(p(xk|xk-1) * Pprev)
        Pint = Pxx * P(ii-1,:)';
        % compute p(zk|xk)
        [Pzx(ii,:) ,pw(ii),qw(ii)]=R_Pzx(Data, Par,xvals,ii);
        % compute P = p(xk|sk) * p(xk|s1,...,sk-1) * p(zk | xk) / p(xk)
        
       if sum(Pzx(ii,:)) ~=0
        p = (Pint'.* Pzx(ii,:)  ) ./ Px;
        else
       p = (Pint'  ) ./ Px;
        end
        
        end
        if sum(p) ~=0
        P(ii,:) =( p/sum(p));
        else
         P(ii,:)=p;
        end
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
        Pxkp1Sk=Pxx * P(ii,:)';
        PSmooth(ii,:)=P(ii,:) .* (Pxx * (PSmooth(ii+1,:)' ./Pxkp1Sk))';
        Pxxkn1YT=  Pxx * P(ii,:)' * ((PSmooth(ii+1,:)' ./ Pxkp1Sk))'  ;
        epj(ii)=xvals * Pxxkn1YT* xvals';
    end
    xkHatS(ii) = xvals*PSmooth(ii,:)';
    xsHatS(ii) = xvals.*xvals*PSmooth(ii,:)' - xkHatS(ii)^2;
end
%% System Identification
options = optimset('Display','off');
%em for pxs
xmu = glmval(Par.XSH.w,  Data.train.SH(1:end-1,:), 'identity')';
Par.XSH.sigma=sqrt(sum( (xvals.^2 * PSmooth(dvd+1:dvd:end,:)') ...
    + xmu.^2 -2.*xmu.*(xvals * PSmooth(dvd+1:dvd:end,:)'))/length(Data.train.t));

fun = @(x)ExpectedQ_PXS_Ns(x,xvals,PSmooth,Data.train.SH,dvd,Par.XSH.sigma);
[bestx,fvalPxs] = fminsearch(fun,Par.XSH.w',options);
if ~isnan(fvalPxs)

Par.XSH.w=bestx';
end
% em for pxz
Par.Z.a0=log(sum(z)/(dt*sum(exp(Par.Z.a2*pw(2:end)'+Par.Z.a3*qw(2:end)').*(exp(xvals )*  PSmooth(2:end,:)'))));

fun = @(x)ExpectedQ_PXZ(x,xkHatS,xvals,PSmooth,pw,qw,zk,dt,Par.Z.a0,Par.Z.a1);
[bestx,fvalPxz] = fminsearch(fun,[Par.Z.a2,Par.Z.a3],options);
if ~isnan(fvalPxz)
% Par.Z.a0=bestx(1);
% Par.Z.a0=bestx(1);
Par.Z.a2=bestx(1);
Par.Z.a3=bestx(2);
end
% 
% em for px0
% options = optimset('Display','final', 'MaxIter',100);
fun = @(x)ExpectedQ_PX0(x,xvals,PSmooth);
[bestx,fvalPx0] = fminsearch(fun,[Par.X0.mu, Par.X0.sigma],options);
if ~isnan(fvalPx0)
Par.X0.mu=bestx(1);
Par.X0.sigma=bestx(2);
end
% 
% % em for pxx
% Par.XX.sigma=sqrt((sum((xvals.^2 * PSmooth(2:end,:)' + (xvals.^2 * PSmooth(1:end-1,:)')- 2*epj(2:end))) - ...
%     sum((xvals.^2 * PSmooth(dvd+1:dvd:end,:)' + (xvals.^2 * PSmooth(1:dvd:end-dvd,:)')- (2*xvals * PSmooth(dvd+1:dvd:end,:)').* (xvals * PSmooth(1:dvd:end-dvd,:)'))))...
%     /(length(Data.train.ti)-(length(Data.train.t))));
% em for pxx
fun = @(x)ExpectedQ_PXX(x,xvals,PSmooth,epj,dvd);
[bestx,fvalPxx] = fminbnd(fun,.2,100);
if ~isnan(fvalPxx)
Par.XX.sigma=bestx(1);
end

Q=-1*sum(fvalPx0+fvalPxz+fvalPxs+fvalPxx);
fprintf ('Ex0=%f, Ezx=%f,EXX=%f, Exs=%f, Q=%f \n',-1*fvalPx0,-1*fvalPxz,-1* fvalPxx,-1*fvalPxs,Q);
fprintf('Sigma XSH=%f,Sigma XX=%f, \n',Par.XSH.sigma,Par.XX.sigma);
if itr/MaxItr-floor(itr/MaxItr) == 0
    tt = Data.train.ti-Data.train.ti(1);
    tx = Data.train.t-Data.train.t(1);
figure
set(gcf,'color','w')

hold on
% h(1) = plot(tx,Data.train.x,'color',[1 0 1],'linewidth',2);
h(1) = plot(tt,xkHat,'color',0*[1 1 1],'linewidth',2);
% he = shadedErrorBar(tt,xkHat,sqrt(xsHat),{'-','color',[0 0 0]},1);
h(2) = plot(tt,xkHatS,'color',[0 1 0],'linewidth',2);
he = shadedErrorBar(tt,xkHatS,sqrt(xsHatS),{'-','color',[0 1 0]},1);
yvalsC = -1*stateBound*ones(size(Data.train.cwords));
h(3) = plot(Data.train.cwords,yvalsC,'or');

yvalsP = 1*stateBound*ones(size(Data.train.pwords));
h(4) = plot(Data.train.pwords,yvalsP,'ob');hold off
hold off
set(gca,'ydir','normal')
ylabel('state','fontsize',22);
xlabel('seconds','fontsize',22);
ylim([-1*stateBound,1*stateBound]);
title(['Estimated State with autoencoder model, h=',num2str(hs)],'fontsize',22);
legend(h,{'kalman-neural','KF-Smooth','Companion','Participant' },'fontsize',14);

pause(1);
end

end
logLikelihoody=nansum( nansum(log( PSmooth(dvd+1:dvd:end,:))));
[bic(bicIndx)]=aicbic(logLikelihoody,length( Par.XSH.w),length(Data.train.x));
end
figure
plot(hs,bic,'ob', 'linewidth',2);
set(gca,'ydir','normal')
ylabel('BIC','fontsize',22);
xlabel('hs','fontsize',22);
%%
save('Apr12','Data','Par','stateBound','xvals','dt','dvd')
run('Reza_DirectDecoder_V1.m')
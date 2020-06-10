%% rebuild neural feature matrix using original data window (train+test)
%Reza Z_k is neureal feature & x_k is position
close all
clear all
clc

load('Data33.mat');
load('penalty.mat');
MSTrain(:,[20,24])=0;

%% preporcessing the data 
gausWindowLength=60;
gausWindowSigma=15;
BestIndxs=R_BestCellsIndx(MSTrain(:,2:63),120)
gausWindow = fspecial('gaussian', gausWindowLength, gausWindowSigma);
gausWindow=mean(gausWindow);
gausWindow=gausWindow/max(gausWindow);
[ySum,yMean, yCells ]= R_smoothNeuralActivity_v2(MSTrain(:,1+BestIndxs), gausWindow);
% yCells=MSTrain(:,1+BestIndxs);
trainPerc=.85;
fprintf('smoothing neural activity with a gaussian window with length= %d and sigma=%2.2f \n', ...
    gausWindowLength, gausWindowSigma);
%% create train dataset
train.data.strInd =  floor(length(MSTrain)*0)+1;
train.data.endInd =  floor(length(MSTrain)*trainPerc);
train.data.zk=yCells(train.data.strInd:train.data.endInd,:);% spiking activity channel
train.data.xk=MSTrain(train.data.strInd:train.data.endInd,64:65);% position channel
train.data.tk=MSTrain(train.data.strInd:train.data.endInd,1); % time channel
fprintf(' creat train data %1.2f-%1.2f \n', 0.0, trainPerc);
%% create test dataset
test.data.strInd =  floor(length(MSTrain)*trainPerc);
test.data.endInd =  floor(length(MSTrain)*1);
test.data.zk=yCells(test.data.strInd:test.data.endInd,:); % spiking activity channel
test.data.xk=MSTrain(test.data.strInd:test.data.endInd,64:65); % position channel
test.data.tk=MSTrain(test.data.strInd:test.data.endInd,1);% time channel
fprintf(' creat test data %1.2f-%1.2f\n', trainPerc,1.0);
%%
v=abs(MSTrain(train.data.strInd:train.data.endInd,67))+abs(MSTrain(train.data.strInd:train.data.endInd,66));
 v=floor(v);
 train.IndxB=find(v >=0 );

train.data.zkv=train.data.zk(train.IndxB,:);
train.data.xkv=train.data.xk(train.IndxB,:);
%%
hy=1:20;
for i=1:length(hy)
fprintf('hy=%d/%d\n',hy(i),max(hy));
train.data.Xdsgny = R_DesignMatrix(train.data.zk,hy(i));
test.data.Xdsgny = R_DesignMatrix(test.data.zk,hy(i));
[by,devy,statsy] = glmfit(train.data.Xdsgny(train.IndxB,:),train.data.xk(train.IndxB,2),'normal');
[yhatt,dylott,dyhitt] = glmval(by,train.data.Xdsgny(train.IndxB,:),'identity',statsy);
logLikelihoody= nansum(log( normpdf( train.data.xk(train.IndxB,2), yhatt, statsy.s)));
[aicy(i),bicy(i)]=aicbic(logLikelihoody,length(by),length(train.data.xk(train.IndxB,:)));
end

hyOpt=hy(find(bicy== min(bicy)));
train.data.Xdsgny = R_DesignMatrix(train.data.zk,hyOpt);
test.data.Xdsgny = R_DesignMatrix(test.data.zk,hyOpt);
[by,devy,statsy] = glmfit(train.data.Xdsgny(train.IndxB,:),train.data.xk(train.IndxB,2),'normal');
[yhatt,dylott,dyhitt] = glmval(by,train.data.Xdsgny(train.IndxB,:),'identity',statsy);
[yhat,dylot,dyhit] = glmval(by,test.data.Xdsgny,'identity',statsy);
%%


hx=1:20;
for i=1:length(hx)
fprintf('hx=%d/%d\n',hx(i),max(hx));
train.data.Xdsgnx = R_DesignMatrix( cat(2,train.data.zk,( 1*yhatt)),hx(i));
test.data.Xdsgnx = R_DesignMatrix(cat(2,test.data.zk,1*(yhat)),hx(i));
[bx,devx,statsx] = glmfit(train.data.Xdsgnx(train.IndxB,:),train.data.xk(train.IndxB,1),'normal');
[xhatt,dxlott,dxhitt] = glmval(bx,train.data.Xdsgnx(train.IndxB,:),'identity',statsx);
logLikelihoodx= nansum(log( normpdf( train.data.xk(train.IndxB,1), xhatt, statsx.s)));
[aicx(i),bicx(i)]=aicbic(logLikelihoodx,length(bx),length(train.data.xk(train.IndxB,:)));
end

hxOpt=hx(find(bicx== min(bicx)));

train.data.Xdsgnx = R_DesignMatrix( cat(2,train.data.zk,( yhatt)),hxOpt);
test.data.Xdsgnx = R_DesignMatrix(cat(2,test.data.zk,(yhat)),hxOpt);
[bx,devx,statsx] = glmfit(train.data.Xdsgnx(train.IndxB,:),train.data.xk(train.IndxB,1),'normal');
[xhat,dxlot,dxhit] = glmval(bx,test.data.Xdsgnx,'identity',statsx);
[xhatt,dxlott,dxhitt] = glmval(bx,train.data.Xdsgnx,'identity',statsx);
%% Plot goodness of fit result
figure
subplot(2,1,1)
plot(hy-1,bicy,'b','linewidth',2);
title('BIC_Y');
xlabel('Hy');

subplot(2,1,2)
plot(hx-1,bicx,'b','linewidth',2);
title('BIC_X');
xlabel('Hx');
%%


figure
hold on
h1 = plot(test.data.xk(:,1),test.data.xk(:,2),'color',[1 0 0]);
h2 = plot(xhat,yhat,'color',[0 1 0],'linewidth',1);
hold off
set(gca,'ydir','normal')
ylabel('state','fontsize',18);
xlabel('seconds','fontsize',18);
title ([' GLM for best Hx=',num2str(hxOpt-1),', Hy=',num2str(hxOpt-1)])
legend([h1,h2],{'true','estimated'},'fontsize',14);
colorbar
fprintf(' fit a GLm model on train data , see the figure \n'); 
%% Boarders Penalty
d_step=5;
x1=floor(mean(test.data.xk(:,1))-4*std(test.data.xk(:,1))):d_step:floor(mean(test.data.xk(:,1))+4*std(test.data.xk(:,1)));
x2=floor(mean(test.data.xk(:,2))-4*std(test.data.xk(:,2))):d_step:floor(mean(test.data.xk(:,2))+4*std(test.data.xk(:,2)));
[X1,X2]  = meshgrid(x1,x2);
for i=1:length(x1)
    for j=1:length(x2)
        Xs(j+(i-1)*length(x2),:) = [x1(i)  x2(j)];
    end
end
% Err=1000*ones(length(Xs),length(train.data.xk));
% for indM=1:length(Xs)
%     for indData=1:length(train.data.xk)
%        Err(indM,indData)= sqrt(distance2(Xs(indM,:),train.data.xk(indData,:)));
%     end
% end

% EErm=min(Err');
% penalty.Indx=(find(EErm >d_step));
    %% build distribution of x empirically
    stdx=std(train.data.xk(4425:4605,1));
    stdy=std(train.data.xk(8280:8555,2));
    Qk=diag([stdx,stdy]);
    Ak = [1 0;0 1];
    Bk = [0.00;0.00];
    Px = mvnpdf([X1(:) X2(:)],mean(train.data.xk),(Qk));
    Px = Px/sum(Px);
    fprintf(' build distribution of x empirically \n'); 
    %%
    nTimeValues = numel(test.data.strInd:test.data.endInd);
    nStateValues = numel(1:length(Px));
    P = nan(nTimeValues,nStateValues);
    xm = nan(nTimeValues,2);
    xss = nan(nTimeValues,2);

    penalty.wtsOut=ones(nStateValues,1)';
    penalty.wtsIn=ones(nStateValues,1)'*1;
    penalty.wtsIn(penalty.Indx)=1;
    penalty.wtsOut(penalty.Indx)=0.000001;
    % compute p(xk|xk-1), i.e. distributions of xk for each possible x value
    % and then normalize columns of p(xk|xk-1)
    Pxx = nan(nStateValues);

    X11=X1(:);
    X22=X2(:);
    %% 
    for ii = 1:nStateValues
        mxx = Ak*[X11(ii),X22(ii)]' + Bk;
        Pxx(ii,:) = mvnpdf(Xs,mxx',Qk/4);
    end
    for ii = 1:nStateValues
        Pxx(:,ii) = Pxx(:,ii)/sum(Pxx(:,ii));
    end
    fprintf(' build xk|xk-1 state empirically \n'); 
    %% Feb6: This is new change
    Pxz_prv = ones(size(Pxx(1,:)))/length(Pxx(1,:)).*penalty.wtsOut;
    xmu=test.data.xk(1,:)';
    hpd=0;
    ss=0;

    for ii = 1:nTimeValues
        fprintf(' calcute D-D output at time index=%d / %d \n',ii,nTimeValues); 
        %% ALI
        %% Feb6: This needs to be changed
        % compute p(xk|hk)~p(xk|z_k-1,hk-1) using variance of x|z (as computed from GLM)
        % this can be done by integral(p(xk|xk-1) * Pxz), where Pxz comes from
        % the previous point
        % what you see here, it is based on 
    %      Px= ones(size(Pxx(1,:)))';
         Px = Pxx * (Pxz_prv'.*penalty.wtsIn');
         Px = Px'/sum(Px');
        if ii == 1
            P(ii,:) = Px .*penalty.wts;
            xm(ii,:) = [X1(:) X2(:)]'*P(ii,:)';
            xss(ii,:) = [X1(:) X2(:)]'.*[X1(:) X2(:)]'*P(ii,:)' - xm(ii)^2;
        else
            % compute p(xk|zk,hk) using variance of x|z (as computed from GLM)
            %%%% need chnge %%%%%%
            [xhn] = glmval(bx,test.data.Xdsgnx(ii,:),'identity',statsx);
            [yhn] = glmval(by,test.data.Xdsgny(ii,:),'identity',statsy);
    %         xmu = [xh,yh];
    %         xmu=test.data.xk(ii,:)+randn([1,2]);
            xh=xhn;
            yh=yhn;
            %%%% %%%%%%%%%%%%%%%
    %         Pxz = mvnpdf([X1(:) X2(:)],xmu,[statsx.s,5;5,statsy.s])';
            Pxz=normpdf(X1(:),xh,(statsy.s))' .* normpdf(X2(:),yh,(statsy.s))';
                
           
            Pxz=Pxz./sum(Pxz);
            Pxz_prv = Pxz; 
            Pxz=Pxz .*penalty.wtsOut;
            % compute p(xk|z1,...,zk-1) = integral(p(xk|xk-1) * Pprev)
            Pint = Pxx * P(ii-1,:)';

            % compute P = p(xk|zk) * p(xk|z1,...,zk-1) / p(xk)

            p = (Pxz) .* (Pint' ./ Px);
%              p = (Pxz) ;
            p=p .*penalty.wtsOut;
            P(ii,:) = p./(sum(p));
    %         

                xm(ii,:) = Xs'*P(ii,:)';
    %                 IndxMax=find(P(ii,:) == max(P(ii,:)));
    %                 xm(ii,:) = (Xs(IndxMax,:));
                xss(ii,:) = Xs'.*Xs'*P(ii,:)' - xm(ii,:)'.^2;
    %         end

        end

        % store x mean/stdev

        %% 
        xa(ii,:)=xmu;

        [se_port,se_ind] = sort(P(ii,:),'descend');
        cse_port         = cumsum(se_port);
        [~,ind] = min(abs(cse_port-0.95));
        dd = sqrt((Xs(se_ind(1:ind),1)-test.data.xk(ii,1)).^2+(Xs(se_ind(1:ind),2)-test.data.xk(ii,2)).^2);
        if min(dd)<d_step
            hpd = hpd +1;
        end
        %%
    end
    hpd = hpd/(test.data.endInd-test.data.strInd+1);
    hpd=hpd*100;
    %%
    %%
    fprintf(' show the D-D model result, see the figure');
    figure
    % set(gcf,'position',[680         558        1117         420])
    set(gcf,'color','w')
    errIndx=1;
    stind=errIndx;
    endind=1537+errIndx;
    Ps=reshape(P,nTimeValues,length(x2),length(x1));
    imagesc(x1,x2, squeeze(sum(Ps(stind:endind,:,:),1)));
    hold on
    h(1) = plot(test.data.xk(stind:endind,1),test.data.xk(stind:endind,2),'color',[1 0 0]);

    h(2) = plot(xm(stind:endind,1),xm(stind:endind,2),'color',0*[1 1 1],'linewidth',2);

    set(gca,'ydir','normal')
    ylabel('Y (cm)','fontsize',22);
    xlabel('X (cm)','fontsize',22);

    legend(h,{'True Traj.','D-D result'},'fontsize',14);
    colorbar
    
    RMSE2d=R_2DRmse(xm(errIndx:end,:),test.data.xk(errIndx:end,:));
    RMSEx=R_1DRmse(xm(errIndx:end,1),test.data.xk(errIndx:end,1));
    RMSEy=R_1DRmse(xm(errIndx:end,2),test.data.xk(errIndx:end,2));
    fprintf(' D-D with Hx=%d,Hy=%d, 2D-RMSE  %1.2f \n',hxOpt,hyOpt,RMSE2d);
    fprintf(' D-D with Hx=%d,Hy=%d, RMSEx %1.2f, RMSEy %1.2f \n',hxOpt,hyOpt,RMSEx,...
       RMSEy);
    fprintf(' D-D with Hx=%d, Hy=%d, HPD  %1.2f \n',hxOpt,hyOpt,hpd);
    %%
    stind=errIndx;
    endind=length(test.data.xk);
    figure
    subplot(211)
    hold on
    imagesc((1:length(stind:endind)).*33,x1,squeeze(sum(Ps(stind:endind,:,:),2))');
    plot((1:length(stind:endind)).*33,test.data.xk(stind:endind,1),'color',1*[0 0 0],'linewidth',4);
    
    plot((1:length(stind:endind)).*33,xm(stind:endind,1),'color',[1 0 0],'linewidth',4)
    hold off
    legend({'True-X','Estimated-X'},'fontsize',14)
    ylabel('X (cm)')
    xlim([1,length(stind:endind).*33])
    ylim([min(x1),max(x1)])
    title(['D-D with Hx=',num2str(hxOpt-1),', Hy=',num2str(hyOpt-1), ', 2D-RMSE', num2str(RMSE2d), ...
        ', RMSEx=', num2str(RMSEx), ', RMSEy=', num2str(RMSEy), ...
        ', HPD=',num2str(hpd)],'fontsize',22)
    subplot(212)
    hold on
    imagesc((1:length(stind:endind)).*33,x2,squeeze(sum(Ps(stind:endind,:,:),3))');

    plot((1:length(stind:endind)).*33,test.data.xk(stind:endind,2),'color',1*[0 0 0],'linewidth',4);
    plot((1:length(stind:endind)).*33,xm(stind:endind,2),'color',1*[1 0 0],'linewidth',4);
    hold off
    legend({'True-Y','Estimated-Y'},'fontsize',14);
    xlabel('Time (ms)','fontsize',14)
    ylabel('Y (cm)','fontsize',14)
    xlim([1,length(stind:endind).*33])
    ylim([min(x2),max(x2)])

    %%
    figure
    subplot(211)
    hold on
%     imagesc(1:length(stind:endind),x1,squeeze(sum(Ps(stind:endind,:,:),2))');

    err = [dxhit(stind:endind)';dxlot(stind:endind)'];
    he = shadedErrorBar((1:length(stind:endind)).*.033,xm(stind:endind,1),sqrt(xss(stind:endind,1)),{'-','color',[0 0 0]},1);
    h1=plot((1:length(stind:endind)).*.033,test.data.xk(stind:endind,1),'color',1*[0 0 0],'linewidth',4);
    h2=plot((1:length(stind:endind)).*.033,xm(stind:endind,1),'color',[1 0 0],'linewidth',4)
    hold off
    legend([h1,h2],{'True-X','Estimated-X'},'fontsize',20)
    ylabel('X (cm)','fontsize',24)
    xlim([1,length(stind:endind).*.033])
    ylim([min(x1),max(x1)])
%     title(['D-D with Hx=',num2str(hxOpt-1),', Hy=',num2str(hyOpt-1), ', 2D-RMSE', num2str(RMSE2d), ...
%         ', RMSEx=', num2str(RMSEx), ', RMSEy=', num2str(RMSEy), ...
%         ', HPD=',num2str(hpd)],'fontsize',22)
    subplot(212)
    hold on
%     imagesc(1:length(stind:endind),x2,squeeze(sum(Ps(stind:endind,:,:),3))');
   
    he = shadedErrorBar((1:length(stind:endind)).*.033,xm(stind:endind,2),sqrt(xss(stind:endind,2)),{'-','color',[0 0 0]},1);
    h1=plot((1:length(stind:endind)).*.033,test.data.xk(stind:endind,2),'color',1*[0 0 0],'linewidth',4);
    h2=plot((1:length(stind:endind)).*.033,xm(stind:endind,2),'color',1*[1 0 0],'linewidth',4);
    hold off
%     legend([h1,h2],{'True-Y','Estimated-Y'},'fontsize',14);
    xlabel('Time (s) ','fontsize',24)
    ylabel('Y (cm)','fontsize',24)
    xlim([1,length(stind:endind).*.033])
    ylim([min(x2),max(x2)])
    %%
    figure

    hold on
    plot(train.data.xk(:,1),train.data.xk(:,2))
    scatter(Xs(penalty.Indx,1),Xs(penalty.Indx,2),'filled')
    hold off
    title('Penalized area','fontsize',22);
    set(gca,'ydir','normal')
    ylabel('Y (cm)','fontsize',22);
    xlabel('X (cm)','fontsize',22);
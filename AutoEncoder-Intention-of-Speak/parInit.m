function Par=parInit(Data,stateBound)
Par.X0.mu=0;
Par.X0.sigma=.5;

Par.XX.Ak = 1;
Par.XX.Bk = 0;
Par.XX.sigma=.5;
Par.XX.mu=0;

Par.Z.a0=0;
Par.Z.a1=1;
Par.Z.a2=-2;
Par.Z.a3=+2;

% [b,dev,stats]=glmfit(Data.train.SH,Data.train.x,'normal');
[M,N]=size(Data.train.SH);
Par.XSH.w=stateBound*ones([N+1,1])/max(max(Data.test.S))/length(Data.bestIdx);
% Par.XSH.w=(b);
Par.XSH.sigma=.5;
Par.XSH.mu=0;


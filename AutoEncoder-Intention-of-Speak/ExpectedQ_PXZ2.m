function UperLimitQ = ExpectedQ_PXZ2(par,xhat,xvals,PXkYT,pw,qw,zk,dt,a0,a1)

%pz parameters
% a0 = par(1);
a2 = par(1);
a3 = par(2);
% a3 = par(4);

dt=1;
Ezx=sum(-1*(dt*exp(a0+a1*xhat(2:end)+a2*pw(2:end)'+a3*qw(2:end)')) + ...
     zk(2:end)'.* (a0+a1*xvals* PXkYT(2:end,:)'+a2*pw(2:end)'+a3*qw(2:end)' +log(dt)));

UperLimitQ=1*(Ezx);
% fprintf ('Ex0=%f, Ezx=%f, Exs=%f, Exx=%f,Exxj=%f, Q=%f \n',Ex0,Ezx,Exs,Exx,Exxj,-1*UperLimitQ);
function UperLimitQ = ExpectedQ_PXZ(par,xhat,xvals,PXkYT,pw,qw,zk,a0,a1)

%pz parameters
% a1 = par(1);
a2 = par(1);
a3 = par(2);
% a3 = par(4);

%     

Ezx=sum(-1*(exp(a0+a1*xhat(2:end)+a2*pw(2:end)'+a3*qw(2:end)')+ ...
        a1*exp(a0+a1*xhat(2:end)+a2*pw(2:end)'+a3*qw(2:end)') .* ((xvals) * PXkYT(2:end,:)')-xhat(2:end)+ ...
    .5*a1^2*((xvals.^2 * PXkYT(2:end,:)') -2*xhat(2:end).*(xvals * PXkYT(2:end,:)')+xhat(2:end).^2) .* exp(a0+a1*xhat(2:end)+a2*pw(2:end)'+a3*qw(2:end)'))+ ...
     zk(2:end)'.* (a0+a1*xvals* PXkYT(2:end,:)'+a2*pw(2:end)'+a3*qw(2:end)' ));
UperLimitQ=-1*(Ezx);
% fprintf ('Ex0=%f, Ezx=%f, Exs=%f, Exx=%f,Exxj=%f, Q=%f \n',Ex0,Ezx,Exs,Exx,Exxj,-1*UperLimitQ);
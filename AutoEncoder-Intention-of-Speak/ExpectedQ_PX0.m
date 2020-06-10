function UperLimitQ = ExpectedQ_PX0(par,xvals,PXkYT)
%Px0 parameters
m0 = par(1);
sigma0 = par(2);


Ex0=sum(-0.5 * log(sigma0^2) -0.5/sigma0^2 * ((xvals.^2 * PXkYT(1,:)') + m0^2 -2*m0*(xvals * PXkYT(1,:)')));

UperLimitQ=-1*(Ex0);
% fprintf ('Ex0=%f, Ezx=%f, Exs=%f, Exx=%f,Exxj=%f, Q=%f \n',Ex0,Ezx,Exs,Exx,Exxj,-1*UperLimitQ);
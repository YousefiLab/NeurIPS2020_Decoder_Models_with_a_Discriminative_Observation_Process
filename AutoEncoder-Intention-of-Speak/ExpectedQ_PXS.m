function UperLimitQ = ExpectedQ_PXS(par,xvals,PXkYT,SH,dvd)
%Pxs parameters
sigmav = par(1);
w=par(2:end);


xmu = glmval(w',  SH(2:end,:), 'identity')';
[M,N]=size((SH(1:end-1,:)));

Exs=sum(-0.5 * log(sigmav^2) -0.5/sigmav^2 * ((xvals.^2 * PXkYT(dvd+1:dvd:end,:)') + ...
    xmu.^2 -2.*xmu.*(xvals * PXkYT(dvd+1:dvd:end,:)')));
UperLimitQ=-1*(Exs);
% fprintf ('Ex0=%f, Ezx=%f, Exs=%f, Exx=%f,Exxj=%f, Q=%f \n',Ex0,Ezx,Exs,Exx,Exxj,-1*UperLimitQ);
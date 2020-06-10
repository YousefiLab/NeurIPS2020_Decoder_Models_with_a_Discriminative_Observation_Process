function UperLimitQ = ExpectedQ_PXS_Ns(par,xvals,PXkYT,SH,dvd,sigmav)
%Pxs parameters

w=par(1:end);


xmu = glmval(w',  SH(1:end-1,:), 'identity')';
[M,N]=size((SH(1:end-1,:)));

Exs=sum(-0.5 * log(sigmav^2) -0.5/(sigmav^2) * ((xvals.^2 * PXkYT(dvd+1:dvd:end,:)') + ...
    xmu.^2 -2.*xmu.*(xvals * PXkYT(dvd+1:dvd:end,:)')));
UperLimitQ=-1*(Exs);
% fprintf ('Ex0=%f, Ezx=%f, Exs=%f, Exx=%f,Exxj=%f, Q=%f \n',Ex0,Ezx,Exs,Exx,Exxj,-1*UperLimitQ);
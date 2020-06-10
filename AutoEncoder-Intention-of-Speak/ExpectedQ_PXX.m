function UperLimitQ = ExpectedQ_PXX(par,xvals,PXkYT,epj,dvd)

% %Pxx parameters
sigmae = par(1);
%pz parameters


Exxj=sum(-0.5 * log(sigmae^2) -0.5/sigmae^2*(xvals.^2 * PXkYT(2:end,:)' + (xvals.^2 * PXkYT(1:end-1,:)')- 2*epj(1:end-1)));
Exx=sum(-0.5 * log(dvd*sigmae^2) -0.5/(dvd*sigmae)^2*(xvals.^2 * PXkYT(dvd+1:dvd:end,:)' + (xvals.^2 * PXkYT(1:dvd:end-dvd,:)')- (2*xvals * PXkYT(dvd+1:dvd:end,:)').* (xvals * PXkYT(1:dvd:end-dvd,:)')));

UperLimitQ=-1*(Exxj-Exx);
% fprintf ('Ex0=%f, Ezx=%f, Exs=%f, Exx=%f,Exxj=%f, Q=%f \n',Ex0,Ezx,Exs,Exx,Exxj,-1*UperLimitQ);
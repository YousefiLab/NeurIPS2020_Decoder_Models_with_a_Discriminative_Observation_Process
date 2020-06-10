function UperLimitQ = ExpectedQ(par,xhat,xvals,PXkYk,PXkYT,Pxxkn1YT,SH,pw,qw)
%Px0 parameters
m0 = par(1);
sigma0 = par(2);
% %Pxx parameters
sigmae = par(3);
%pz parameters
a0 = par(4);
a1 = par(5);
a2 = par(6);
a3 = par(7);
%Pxs parameters
sigmav = par(8);
w=par(9:end);
zk=zeros(size(pw));
zk=pw+qw;
zk(find(zk >0))=1;
Ex0=sum(-0.5 * log(sigma0^2) -0.5/sigma0^2 * ((xvals.^2 * PXkYT(1,:)') + m0^2 -2*m0*(xvals * PXkYT(1,:)')));

Exx=sum(-0.5 * log(sigmae^2) -0.5/sigmae^2*(xvals.^2 * PXkYT(2:end,:)' + (xvals.^2 * PXkYT(1:end-1,:)')- (2*xvals * PXkYT(2:end,:)').* (xvals * PXkYT(1:end-1,:)')));
Exxj=sum(-0.5 * log(sigmae^2) -0.5/sigmae^2*(xvals.^2 * Pxxkn1YT(2:end,:)' + (xvals.^2 * Pxxkn1YT(1:end-1,:)')- (2*xvals * Pxxkn1YT(2:end,:)').* (xvals * Pxxkn1YT(1:end-1,:)')));

Ezx=sum(-exp(a0+a1*xhat(2:end)+a2*pw(2:end)'+a3*qw(2:end)')+ ...
    .5*a1^2*((xvals.^2 * PXkYT(2:end,:)') -2*xhat(2:end).*(xvals * PXkYT(2:end,:)')+xhat(2:end).^2) .* exp(a0+a1*xhat(2:end)+a2*pw(2:end)'+a3*qw(2:end)') ...
    + zk(2:end)'.* (a0+a1*xvals* PXkYT(2:end,:)'+a2*pw(2:end)'+a3*qw(2:end)'));

xmu = glmval(w',  SH(2:end,:), 'identity')';
[M,N]=size((SH(2:end,:)));
xmu=upsample(xmu,floor(length(xhat)/M));
Exs=sum(-0.5 * log(sigmav^2) -0.5/sigmav^2 * ((xvals.^2 * PXkYT(2:end,:)') + ...
    xmu.^2 -2.*xmu.*(xvals * PXkYT(2:end,:)')));
UperLimitQ=-1*(Ex0+Exs+Ezx+Exxj-Exx);
% fprintf ('Ex0=%f, Ezx=%f, Exs=%f, Exx=%f,Exxj=%f, Q=%f \n',Ex0,Ezx,Exs,Exx,Exxj,-1*UperLimitQ);
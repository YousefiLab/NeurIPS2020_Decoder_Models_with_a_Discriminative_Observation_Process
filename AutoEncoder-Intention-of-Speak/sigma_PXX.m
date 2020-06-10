function sigma = sigma_PXX(xvals,PXkYT,epj,dvd)




Exxj=sum((xvals.^2 * PXkYT(2:end,:)' + (xvals.^2 * PXkYT(1:end-1,:)')- 2*epj(2:end)));
Exx=sum((xvals.^2 * PXkYT(dvd+1:dvd:end,:)' + (xvals.^2 * PXkYT(dvd:dvd:end-1,:)')- (2*xvals * PXkYT(dvd+1:dvd:end,:)').* (xvals * PXkYT(dvd:dvd:end-1,:)')));

sigma= sqrt(Exxj-Exx)
% fprintf ('Ex0=%f, Ezx=%f, Exs=%f, Exx=%f,Exxj=%f, Q=%f \n',Ex0,Ezx,Exs,Exx,Exxj,-1*UperLimitQ);
function I=R_BestCellsIndx(y,thr)

I=find(sum(y) >= thr);

function rmse2d=R_2DRmse(x,xhat)
a=zeros(length(x),1);
for i =1: length(x)
    
    a(i)=distance2(x(i,:), xhat(i,:));
end
rmse2d=sqrt(mean(a));

function rmse2d=R_2DRmseMin(x,xhat)
a=zeros(length(x),1);
for i =1: length(x)
    
    a(i)=sqrt(distance2(x(i,:), xhat(i,:)));
end
rmse2d=(min(a));

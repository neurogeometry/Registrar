% this function finds the smalest pad=N-N0 which will factorize N 

function pad=factorize(N0)

factors=[2,3,5,7,11,13];
max_powers=ceil(log(N0)./log(factors));

[X1,X2,X3,X4,X5,X6]=ndgrid(0:max_powers(1),0:max_powers(2),0:max_powers(3),0:max_powers(4),0:max_powers(5),0:max_powers(6));
temp=factors(1).^X1.*factors(2).^X2.*factors(3).^X3.*factors(4).^X4.*factors(5).^X5.*factors(6).^X6-N0;
temp(temp<0)=inf;
[pad,~]=min(temp(:));
%[n1,n2,n3,n4,n5]=ind2sub(size(temp),ind);
%n1=n1-1; n2=n2-1; n3=n3-1; n4=n4-1;; n5=n5-1;
%[n1,n2,n3,n4,n5]

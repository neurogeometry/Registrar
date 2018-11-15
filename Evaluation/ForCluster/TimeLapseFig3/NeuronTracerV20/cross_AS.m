% this function calculates the cross product
% a, b, and c are in a Nx3 format

function c=cross_AS(a,b)

c=zeros(size(a));
c(:,1)=a(:,2).*b(:,3)-b(:,2).*a(:,3);
c(:,2)=-a(:,1).*b(:,3)+b(:,1).*a(:,3);
c(:,3)=a(:,1).*b(:,2)-b(:,1).*a(:,2);
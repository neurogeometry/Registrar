% This function converts subscript to index format

function ind =sub2ind_AS(sizeIm,x,y,z)

if max(x)>sizeIm(1) | max(y)>sizeIm(2)%| max(z)>sizeIm(3)
    error('Out of range subscript')
else
    ind=x+(y-1).*sizeIm(1)+(z-1).*(sizeIm(1)*sizeIm(2));
end

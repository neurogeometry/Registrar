% This function converts index to subscript format

function [x,y,z] =ind2sub_AS(sizeIm,ind)

if max(ind)>prod(sizeIm) 
    error('Out of range index')
else
    z=ceil(ind./(sizeIm(1)*sizeIm(2)));
    y=ceil((ind-(z-1).*(sizeIm(1)*sizeIm(2)))./sizeIm(1));
    x=ind-(z-1).*(sizeIm(1)*sizeIm(2))-(y-1).*sizeIm(1); 
end
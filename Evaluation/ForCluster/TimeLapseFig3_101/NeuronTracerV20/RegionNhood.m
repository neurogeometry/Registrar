% This function finds unique indices of all neighbors of sub_ind and 
% returns them in a 1D array. 
% sub_ind can be in index or subscript format
% ii, jj, kk represent the structural element.

function NHOOD = RegionNhood(sub_ind,sizeIm,ii,jj,kk)

if size(sub_ind,2)==1
    [x,y,z]=ind2sub_AS(sizeIm,sub_ind);
elseif size(sub_ind,2)==3
    x=sub_ind(:,1); 
    y=sub_ind(:,2);  
    z=sub_ind(:,3); 
end

X=ii*ones(1,length(x))+ones(length(ii),1)*x';
Y=jj*ones(1,length(y))+ones(length(jj),1)*y';
Z=kk*ones(1,length(z))+ones(length(kk),1)*z';

temp=(X(:)>=1 & X(:)<=sizeIm(1) & Y(:)>=1 & Y(:)<=sizeIm(2) & Z(:)>=1 & Z(:)<=sizeIm(3));

NHOOD = X(temp)+(Y(temp)-1).*sizeIm(1)+(Z(temp)-1).*(sizeIm(1)*sizeIm(2));
NHOOD = NHOOD(:);
NHOOD(NHOOD<1 | NHOOD>prod(sizeIm))=[];
NHOOD = unique(NHOOD);

% This function returns the X1, X2, Y1, Y2, Z1, Z2, and Color (RGB) arrays for tree 
% structures contained in AM. 
% The function works with labeled or not labeled AM.
% AM can be directed or undirected. 
% The labels don't have to be consecutive.

function [X1,X2,Y1,Y2,Z1,Z2,C]=PlotAM_XYZC(AM,r)

r=r-0.5; % r ranges from 0.5 to sizeIm+0.5 in Matlab and [0 sizeIm] in Java and SWC

AM = max(AM,AM');
AM = triu(AM);

Labels=unique(full(AM(AM(:)>0)));
L=length(Labels);

if L==1
    C_tree=[1,0,0]; %'m';
else   
    C_tree=hsv(3*L);
    C_tree=C_tree(randperm(3*L),:);
    
    ind=(sum((C_tree-ones(size(C_tree,1),1)*[1,0,0]).^2,2).^0.5>0.3);
    C_tree=C_tree(ind,:);
    C_tree=C_tree(1:L,:);
    C_tree(1,:)=[1,0,0];
end

X1=zeros(nnz(AM),1); X2=X1; Y1=X1; Y2=X1; Z1=X1; Z2=X1;
C=zeros(nnz(AM),3);
count=0;
for f=1:L
    [i,j]=find(AM==Labels(f));
    Y1(count+1:count+length(i))=r(i,1);
    Y2(count+1:count+length(i))=r(j,1);
    X1(count+1:count+length(i))=r(i,2);
    X2(count+1:count+length(i))=r(j,2);
    Z1(count+1:count+length(i))=r(i,3);
    Z2(count+1:count+length(i))=r(j,3);
    C(count+1:count+length(i),:)=ones(length(i),1)*C_tree(f,:);
    count=count+length(i);
end
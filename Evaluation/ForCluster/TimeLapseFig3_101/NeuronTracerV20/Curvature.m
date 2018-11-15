% This function calculates the curvature for every tree contained in AMlbl
% flag=1 corresponds to the geometric calculation of curvature, Ki=sin/[(l1+l2)/2] 
% flag=2 is the calculation based on |r'*r''|/|r'|^3, Ki=l1*l2*sin/[(l1+l2)/2]^3 

function [Kmax,Ktotal,Kmean,Kcv]=Curvature(AMlbl,r,flag)

Labels=unique(AMlbl(AMlbl>0));
Kmax=zeros(1,max(Labels));
Ktotal=zeros(1,max(Labels));
Kmean=zeros(1,max(Labels));
Kcv=zeros(1,max(Labels));

for i=1:length(Labels)
    [e1,e2]=find(AMlbl==Labels(i));
    e12=unique([e1;e2]);
    AM_tree=AMlbl(e12,e12);
    AM_tree(AM_tree~=Labels(i))=0;
    AM_tree=spones(AM_tree+AM_tree');
    r_tree=r(e12,:);
    
    l1=nan(size(AM_tree,1),1);
    l2=nan(size(AM_tree,1),1);
    sin_theta=nan(size(AM_tree,1),1);
    for k=1:size(r_tree,1)
        temp=find(AM_tree(k,:));
        if length(temp)==2
            temp1=r_tree(temp(2),:)-r_tree(k,:);
            temp2=r_tree(k,:)-r_tree(temp(1),:);
            l1(k)=sum(temp1.^2)^0.5;
            l2(k)=sum(temp2.^2)^0.5;
            sin_theta(k)=real((1-(sum(temp1.*temp2)/l1(k)/l2(k))^2)^0.5);
        end
    end
    
    if flag==1
        Kmax(Labels(i))=max(sin_theta./((l1+l2)./2));
        Ktotal(Labels(i))=nansum(sin_theta);
        Kmean(Labels(i))=nanmean(sin_theta./((l1+l2)./2));
        Kcv(Labels(i))=nanstd(sin_theta./((l1+l2)./2))/Kmean(Labels(i));
    elseif flag==2
        Kmax(Labels(i))=max(l1.*l2.*sin_theta./((l1+l2)./2).^3);
        Ktotal(Labels(i))=nansum(l1.*l2.*sin_theta./((l1+l2)./2).^2);
        Kmean(Labels(i))=nanmean(l1.*l2.*sin_theta./((l1+l2)./2).^3);
        Kcv(Labels(i))=nanstd(l1.*l2.*sin_theta./((l1+l2)./2).^3)/Kmean(Labels(i));
    else
        error('incorrect flag')
    end
end


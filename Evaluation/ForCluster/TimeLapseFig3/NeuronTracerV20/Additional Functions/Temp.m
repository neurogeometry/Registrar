function Untitled

a=[5.*ones(1,1),0.*ones(1,5),0.*ones(1,5)];
a=a+1.*rand(1,length(a));


[~,V_ind]=find_best_cut(a);
%sa=sort(a);
%sV_ind=find_best_cut(sa);

figure, hold on
plot(a,'b*-')
plot([V_ind,V_ind]+0.5,[0,3],'g-')
%plot([1,length(sa)],(sa(sV_ind)+sa(sV_ind+1))/2.*[1,1],'r-')

% [V_min_l,V_ind_l]=find_best_cut(a(1:V_ind))
% [V_min_r,V_ind_r]=find_best_cut(a(V_ind+1:end))
% V_ind_r=V_ind_r+V_ind;
% 
% plot([V_ind_l,V_ind_l]+0.5,[0,3],'r-')
% plot([V_ind_r,V_ind_r]+0.5,[0,3],'r-')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [minC,ind]=find_best_cut(a)
V=zeros(1,length(a));
C=zeros(1,length(a));
for i=1:length(a)-1
    V(i)=(i*var(a(1:i),1)+(length(a)-i)*var(a(i+1:end),1))/length(a);
    temp=(mean(a(1:i))-mean(a(i+1:end)));
    temp(temp<0)=0;
    C(i)=temp;
end
V(length(a))=var(a,1);
C(length(a))=mean(a);
Cost=V.^0.5./C./(1:length(a)).^0.5;
[minC,ind]=min(Cost);


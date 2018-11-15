% This function calculates the distance from any point in AM to the 
% nearest branch or end point. 

function D = Dist_to_bpendp(AM,r)

[it,jt]=find(AM);
N=length(AM);
lll=AM; lll(AM>0)=sum((r(it,:)-r(jt,:)).^2,2).^0.5;
intermp=find(sum(AM)==2);

D=zeros(N,1);
D(intermp)=inf;
for i=1:length(intermp)
    AMtemp=AM;
    ind12=find(AM(intermp(i),:));
    ind1=ind12(1);
    ind2=ind12(2);
    d1=lll(ind1,intermp(i));
    d2=lll(ind2,intermp(i));
    AMtemp(intermp(i),ind12)=0;
    AMtemp(ind12,intermp(i))=0;
    count=1;
    while sum(AMtemp(:,ind1))==1 && sum(AMtemp(:,ind2))==1 && count<=5
        count=count+1;
        ind1new=find(AMtemp(ind1,:));
        ind2new=find(AMtemp(ind2,:));
        AMtemp(ind1,ind1new)=0;
        AMtemp(ind1new,ind1)=0;
        AMtemp(ind2,ind2new)=0;
        AMtemp(ind2new,ind2)=0;
        d1=d1+lll(ind1,ind1new);
        d2=d2+lll(ind2,ind2new);
        ind1=ind1new;
        ind2=ind2new;
    end
    if sum(AMtemp(:,ind1))~=1 && sum(AMtemp(:,ind2))~=1
        D(intermp(i))=min(d1,d2);
    elseif sum(AMtemp(:,ind1))~=1
        D(intermp(i))=d1;
    elseif sum(AMtemp(:,ind1))~=1
        D(intermp(i))=d2;
    end
end



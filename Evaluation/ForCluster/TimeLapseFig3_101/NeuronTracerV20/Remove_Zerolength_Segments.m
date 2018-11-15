% This function eliminates zero length segments from AM   

function [AM r] = Remove_Zerolength_Segments(AM,r)

AM=spones(triu(AM));
AMsym=(AM+AM');
[i j] = find(AM);
d = sum((r(j,:)-r(i,:)).^2,2).^0.5; 

ShortSegs = (d==0);
pre=i(ShortSegs);
post=j(ShortSegs);

if ~isempty(pre)
    rem_vert=unique([pre;post]);
    rem_vert((sum(AMsym(:,rem_vert)))~=2)=[];
    rem_vert(sum(AM(rem_vert,rem_vert))>0)=[];
    
    temp=find(AMsym(:,rem_vert));
    [prepost,~]=ind2sub(size(AMsym(:,rem_vert)),temp);
    pre=prepost(1:2:end-1);
    post=prepost(2:2:end);
    ind=sub2ind(size(AM),pre,post);
    AM(ind)=1;
    
    AM(rem_vert,:)=[];
    AM(:,rem_vert)=[];
    
    r(rem_vert,:)=[];
end
AM=AM+AM';
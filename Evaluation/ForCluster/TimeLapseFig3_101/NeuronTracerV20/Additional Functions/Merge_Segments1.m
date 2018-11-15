% This function merges all segments shorter than 1/2/ppm with one of their 
% neighboring segments. ppm is the number of points per micrometer.

function [AM r R] = Merge_Segments1(AM,r,R,ppm)

AM=spones(triu(AM));
AMsym=(AM+AM');
type=sum(AMsym,2);
[pre post] = find(AM);
d = sum((r(pre,:)-r(post,:)).^2,2).^0.5; 

ShortSegs = (d<=1/ppm/2 & type(pre)==2 & type(post)==2);
pre=pre(ShortSegs);
post=post(ShortSegs);

if ~isempty(pre)
    temp=[pre,post];
    for i=2:length(pre)
        if nnz(temp(1:i-1,:)==pre(i) | temp(1:i-1,:)==post(i))>0
            temp(i,:)=[0,0];
        end
    end
    
    rem_segs=(sum(temp,2)==0);
    pre(rem_segs)=[];
    post(rem_segs)=[];
    prepost_ind=sub2ind(size(AM),pre,post);
    AM(prepost_ind)=0;
    
    [pre1,post1]=find(AM(post,:));
    [post2,pre2]=find(AM(:,post));
    ind1=sub2ind(size(AM),pre(pre1),post1);
    ind2=sub2ind(size(AM),pre(pre2),post2);
    AM([ind1;ind2])=1;
    AM(post,:)=[];
    AM(:,post)=[];
    
    r(pre,:)=(r(pre,:)+r(post,:))./2;
    R(pre)=(R(pre)+R(post))./2;
    r(post,:)=[];
    R(post)=[];
end
AM=AM+AM';

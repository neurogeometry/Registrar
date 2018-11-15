% This function divides all segments longer than 2/ppm in two.
% ppm is the number of points per micrometer

function [AM r R] = Divide_Segments(AM,r,R,ppm)

AM=spones(triu(AM));
[i j] = find(AM);
d = sum((r(j,:)-r(i,:)).^2,2).^0.5; 

LongSegs = (d>=2/ppm);
pre=i(LongSegs);
post=j(LongSegs);

if ~isempty(pre)
    L=length(AM);
    L_add=length(pre);
    
    ind=sub2ind([L,L],pre,post);
    AM(ind)=0;
    AM(L+L_add,L+L_add)=0;
    ind=sub2ind([L+L_add,L+L_add],pre,((L+1):(L+L_add))');
    AM(ind)=1;
    ind=sub2ind([L+L_add,L+L_add],post,((L+1):(L+L_add))');
    AM(ind)=1;
    
    r=[r;(r(pre,:)+r(post,:))./2];
    R=[R;(R(pre)+R(post))./2];
end
AM=AM+AM';





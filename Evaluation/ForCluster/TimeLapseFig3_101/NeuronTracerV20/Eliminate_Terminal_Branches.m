% This function eliminates all fitst-order terminal branches that are shorter than L_thr.
% The input and output are undirected.
% If flag1=1, the final matrix is reduced. 
% If flag2=1, the ROOT is not eliminated.

function [AMlbl_short r_short BranchN] = Eliminate_Terminal_Branches(AMlbl,r,L_thr,flag1,flag2)

r_short=r;
%R_short=R;
AMlbl_short=AMlbl;
AM=spones(AMlbl);
endp=(sum(AM)==1);
endp_Label=sum(AMlbl(:,endp));

BranchLengths=BranchLengthsAM(AMlbl,r);

L=endp_Label(BranchLengths(endp_Label)<L_thr);
if flag2==1
    L(L==1)=[];
end

for i=1:length(L)
    AMlbl_short(AMlbl==L(i))=0;
end
BranchN=length(unique(AMlbl(AMlbl>0)))-length(L);

% cut and relabel AMlbl_short, r_short, and R_short
cut=(sum(AMlbl_short)==0);
if nnz(cut)>0 && flag1==1
    AMlbl_short(cut,:)=[];
    AMlbl_short(:,cut)=[];
    r_short(cut,:)=[];
    %R_short(cut)=[];
    AMlbl_short = LabelBranchesAM(AMlbl_short>0);
end
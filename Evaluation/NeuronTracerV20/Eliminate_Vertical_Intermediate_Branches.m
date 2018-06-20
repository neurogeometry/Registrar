% This function finds and eliminates short vertical intermediate branches. 

function [AMlbl_short,r_short] = Eliminate_Vertical_Intermediate_Branches(AMlbl,r,L_thr)

AM=spones(AMlbl);

% find short intermediate branches 
BranchLengths=BranchLengthsAM(AMlbl,r);
ReduceLabels=find(BranchLengths<L_thr);
ind=zeros(1,length(ReduceLabels));
e1=zeros(1,length(ReduceLabels));
e2=zeros(1,length(ReduceLabels));
for i=1:length(ReduceLabels)
    [e1_temp ~]=find(AMlbl==ReduceLabels(i));
    e1_temp=e1_temp(sum(AM(:,e1_temp))>=3);
    if length(e1_temp)==2
        e1(i)=e1_temp(1);
        e2(i)=e1_temp(2);
        r_temp=r(e1(i),:)-r(e2(i),:);
        if r_temp(3)/(r_temp(1)^2+r_temp(2)^2)^0.5>2
            ind(i)=1;
        end
    end
end
ReduceLabels=ReduceLabels(ind>0);
e1=e1(ind>0);
e2=e2(ind>0);

AMlbl_short=AMlbl;
r_short=r;
%R_short=R;
for i=1:length(ReduceLabels)  
    AMlbl_short(AMlbl_short==ReduceLabels(i))=0;    
    AMlbl_short(e1(i),:)=0; AMlbl_short(:,e1(i))=0;
    AMlbl_short(e2(i),:)=0; AMlbl_short(:,e2(i))=0;
end

% cut and relabel AMlbl_short, r_short and R_short
cut=(sum(AMlbl_short)==0);
if nnz(cut)>0
    AMlbl_short(cut,:)=[];
    AMlbl_short(:,cut)=[];
    r_short(cut,:)=[];
    %R_short(cut)=[];
    AMlbl_short = LabelBranchesAM(AMlbl_short>0);
end



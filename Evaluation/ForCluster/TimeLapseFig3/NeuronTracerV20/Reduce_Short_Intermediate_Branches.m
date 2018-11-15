% This function finds and reduces short intermediate branches. Vertex pairs for 
% every short intermediate branches are replaced with a single vertex
% located at the CM of that branch.

function [AMlbl_short,r_short] = Reduce_Short_Intermediate_Branches(AMlbl,r,L_thr)

AM=spones(AMlbl);

% find short intermediate branches 
BranchLengths=BranchLengthsAM(AMlbl,r);
ReduceLabels=find(BranchLengths<L_thr);
ind=zeros(1,length(ReduceLabels));
e1=zeros(1,length(ReduceLabels));
e2=zeros(1,length(ReduceLabels));
for i=1:length(ReduceLabels)
    [e1_temp e2_temp]=find(AMlbl==ReduceLabels(i));
    e1_temp=e1_temp(sum(AM(:,e1_temp))>=3);
    if length(e1_temp)==2
        e1(i)=e1_temp(1);
        e2(i)=e1_temp(2);
        ind(i)=1;
    end
end
ReduceLabels=ReduceLabels(ind>0);
BranchLengths=BranchLengths(ReduceLabels);
e1=e1(ind>0);
e2=e2(ind>0);

% eliminate contacting intermediate branches
temp=sort([e1,e2]);
rem_ind=[];
temp_ind=find(temp(1:end-1)==temp(2:end));
if ~isempty(temp_ind)
    V=temp(temp_ind);
    for i=1:length(V)
        temp_ind=find([e1,e2]==V(i));
        temp_ind(temp_ind>length(e1))=temp_ind(temp_ind>length(e1))-length(e1);
        [temp,temp_ind1]=sort(BranchLengths(temp_ind));
        rem_ind=[rem_ind,temp_ind(temp_ind1(2))];
    end
end
ReduceLabels(rem_ind)=[];
e1(rem_ind)=[];
e2(rem_ind)=[];

AMlbl_short=AMlbl;
r_short=r;
%R_short=R;
for i=1:length(ReduceLabels)  
    AMlbl_short(AMlbl_short==ReduceLabels(i))=0;
    
    % reconnect
    ends1=find(AMlbl_short(e1(i),:));
    ends2=find(AMlbl_short(e2(i),:));
    AMlbl_short(e1(i),:)=0; AMlbl_short(:,e1(i))=0;
    AMlbl_short(e2(i),:)=0; AMlbl_short(:,e2(i))=0;
    AMlbl_short(ends1,end+1)=1; AMlbl_short(end+1,ends1)=1;
    AMlbl_short(ends2,end)=1; AMlbl_short(end,ends2)=0;
    r_short = [r_short; (r(e1(i),:)+r(e2(i),:))./2];
    %R_short = [R_short; (R(e1(i))+R(e2(i)))./2];
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



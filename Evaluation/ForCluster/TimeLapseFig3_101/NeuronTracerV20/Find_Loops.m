% This function finds the loop verticies in AM.
% Loop_vertices contains the vertices of all branches involved in loops.
% There should be no loops formed by adjacent vertices in AM.
% AM is symmetrized and labeld
% Numbers of loop clusters and simple loops are displayed

function Loop_vertices = Find_Loops(AMlbl)

% label AM
AMlbl=LabelBranchesAM(AMlbl);
AM=spones(AMlbl);
L=unique(AMlbl(AMlbl>0));
Loop_labels=[];

% (1) Reduce the structure by eliminating intermediate vertices
% Terminal and isolated loops are accounted for
AMlbl_short = sparse(size(AMlbl,1),size(AMlbl,2));
Isolated_Loops=0;
Terminal_Loops=0;
for i=1:length(L)
    [e1 ~]=find(AMlbl==L(i));
    e1=unique(e1);    
    temp=e1(sum(AM(:,e1))==1 | sum(AM(:,e1))>=3);
    
    if isempty(temp) % isolated loop
        Loop_labels=[Loop_labels;full(L(i))];
        Isolated_Loops=Isolated_Loops+1;
    elseif length(temp)==1 % terminal loop
        Loop_labels=[Loop_labels;full(L(i))];
        Terminal_Loops=Terminal_Loops+1;
    else
        if AMlbl_short(temp(1),temp(2))==0
            AMlbl_short(temp(1),temp(2))=L(i);
        elseif AMlbl_short(temp(2),temp(1))==0
            AMlbl_short(temp(2),temp(1))=L(i);
        else
            Loop_labels=[Loop_labels;full(L(i))];
        end
    end
end

% (2) Cut AMlbl_short
cut=(sum(AMlbl_short)==0 & sum(AMlbl_short,2)'==0);
if nnz(cut)>0
    AMlbl_short(cut,:)=[];
    AMlbl_short(:,cut)=[];
end

% (3) Remove terminal branches not containing loops
BranchN=-1;
BranchN_new=length(unique(AMlbl_short(AMlbl_short>0)));
while BranchN_new~=BranchN && BranchN_new>0
    BranchN=BranchN_new;
    AM=AMlbl_short>0;
    endp=(sum(AM>0,1)+sum(AM>0,2)'==1);
%     AM=(AMlbl_short | AMlbl_short');
%     endp=(sum(AM,1)==1);
    endp_Label=sum(AMlbl_short(:,endp),1)+sum(AMlbl_short(endp,:),2)';
    
    for i=1:length(endp_Label)
        AMlbl_short(AMlbl_short==endp_Label(i))=0;
    end
    BranchN_new=length(unique(AMlbl_short(AMlbl_short>0)));
end

% (4) Cut AMlbl_short
cut=(sum(AMlbl_short)==0 & sum(AMlbl_short,2)'==0);
if nnz(cut)>0
    AMlbl_short(cut,:)=[];
    AMlbl_short(:,cut)=[];
end

% (5) Cut intermediate branches connnecting loop clusters
L=unique(AMlbl_short(AMlbl_short>0));
[ClusterN,~]=FindClustersAM(AMlbl_short);
for i=1:length(L)
    AMtemp=AMlbl_short;
    AMtemp(AMlbl_short==L(i))=0;
    stay=(sum(AMtemp,1) | sum(AMtemp,2)');
    AMtemp=AMtemp(stay,stay);
    [ClusterNtemp,~]=FindClustersAM(AMtemp);
    if ClusterNtemp==ClusterN
        Loop_labels=[Loop_labels;full(L(i))];
    end
end

% (6) Find loop verticies
Loop_vertices=[];
AMlbl_loops=sparse(size(AMlbl_short,1),size(AMlbl_short,2));
for i=1:length(Loop_labels)
    [e1 ~]=find(AMlbl==Loop_labels(i));
    Loop_vertices=[Loop_vertices;unique(e1)];
    
    AMlbl_loops=AMlbl_loops+(AMlbl_short==Loop_labels(i)).*Loop_labels(i);
end
Loop_vertices=unique(Loop_vertices);
cut=(sum(AMlbl_loops)==0 & sum(AMlbl_loops,2)'==0);
if nnz(cut)>0
    AMlbl_loops(cut,:)=[];
    AMlbl_loops(:,cut)=[];
end

% Number of individual loops = ClusterNumber + Edges - Vertices;
ClusterNumber=FindClustersAM(AMlbl_loops);
LoopNumber = ClusterNumber+length(Loop_labels)-length(AMlbl_loops);
display(['There are ',num2str(ClusterNumber+Isolated_Loops+Terminal_Loops),' loop clusters containing ', num2str(LoopNumber),' simple loops.'])




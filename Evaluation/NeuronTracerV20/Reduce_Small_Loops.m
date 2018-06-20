% This function finds and eliminates small loops from AMlbl. Every loop is
% replaced with a new vertex located at the loop's center of intesnsity.
% Branches longer than MinLoopSize are cut to find small loops.
% ClusterNumber is the number of loop clusters, LoopNumber is the total number of
% loops.

function [AMlblout rout ClusterNumber LoopNumber] = Reduce_Small_Loops(AMlbl,r,MinLoopSize)

% Cut long branches to eliminate large loops
BranchLengths=BranchLengthsAM(AMlbl,r);
CutLabels=find(BranchLengths>MinLoopSize);
AMlbl_short=AMlbl;
AMlbl_topology = Eliminate_Intermediate_Vertices(AMlbl_short); % reduced matrix
r_short=r;
%R_short=R;
for i=1:length(CutLabels)
    AMlbl_short(AMlbl_short==CutLabels(i))=0;
    AMlbl_topology(AMlbl_topology==CutLabels(i))=0;
end

% Cut AMlbl_topology
cut=(sum(AMlbl_topology)==0 & sum(AMlbl_topology,2)'==0);
if nnz(cut)>0
    AMlbl_topology(cut,:)=[];
    AMlbl_topology(:,cut)=[];
end

% Remove terminal branch structures not containing loops
BranchN_new=length(unique(AMlbl_short(AMlbl_short>0)));
BranchN=0;
while BranchN_new~=BranchN
    BranchN=BranchN_new;
    [AMlbl_short r_short BranchN_new] = Eliminate_Terminal_Branches(AMlbl_short,r_short,Inf,0,0); 
end

% Cut intermediate branches connnecting loop clusters
AMlbl_loops=AMlbl_short; % This matrix only contains isolated loop clusters 
L=unique(AMlbl_topology(AMlbl_topology>0));
[ClusterN,~]=FindClustersAM(AMlbl_topology);
for i=1:length(L)
    AMtemp=AMlbl_topology;
    AMtemp(AMtemp==L(i))=0;
    stay=(sum(AMtemp,1) | sum(AMtemp,2)');
    AMtemp=AMtemp(stay,stay);
    [ClusterNtemp,~]=FindClustersAM(AMtemp);
    if ClusterNtemp>ClusterN
        AMlbl_loops(AMlbl_loops==L(i))=0;
    end
end

[ClusterNumber,ClusterVertices] = FindClustersAM(AMlbl_loops); 
% Number of individual loops = ClusterNumber + Edges - Vertices; 
LoopNumber = ClusterNumber + nnz(AMlbl_loops)/2 - nnz(sum(AMlbl_loops,1)); 

% Replace individual clusters with single CM points
AMlblout = AMlbl; 
rout=r;
%Rout=R;

for i=1:ClusterNumber
    ClusterAM=AMlbl(ClusterVertices{i},ClusterVertices{i});
    ClusterLabels=unique(ClusterAM(ClusterAM>0));
    for j=1:length(ClusterLabels)
        AMlblout(AMlblout==ClusterLabels(j))=0;
    end
    
    % Find extenal vertices
    Vtemp=ClusterVertices{i}(sum(AMlblout(:,ClusterVertices{i}))>0);
    
    % Replace all Cluster Vertices with CM Verex
    rout = [rout; mean(r(Vtemp,:),1)];
    %Rout = [Rout; mean(R(Vtemp))];
    
    AMlblout(end+1,end+1)=0; 
    AMlblout(Vtemp,size(AMlblout,2)) = sum(AMlblout(:,Vtemp));
    AMlblout(size(AMlblout,2),Vtemp) = sum(AMlblout(:,Vtemp));
end

% cut and relabel AMlblout, rout, and Rout
cut=(sum(AMlblout)==0);
if nnz(cut)>0
    AMlblout(cut,:)=[];
    AMlblout(:,cut)=[];
    rout(cut,:)=[];
    %Rout(cut)=[];
    AMlblout = LabelBranchesAM(AMlblout>0);
end

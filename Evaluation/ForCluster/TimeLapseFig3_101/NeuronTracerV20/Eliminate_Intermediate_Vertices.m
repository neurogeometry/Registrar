% This function eliminates all intermediate vertices from undirected AMlbl and returns
% the result in a labeled directed matrix of the same size. The output has to
% be directed to encode two vertex loops.

function AMlbl_topology = Eliminate_Intermediate_Vertices(AMlbl)

AMlbl_topology = sparse(size(AMlbl,1),size(AMlbl,2));
AM=spones(AMlbl);
L=unique(AMlbl(AM>0));

% connect branch or end vertices directionally using the branch label
for i=1:length(L)
    [e1 ~]=find(AMlbl==L(i));
    temp=e1(sum(AM(:,e1))==1 | sum(AM(:,e1))>=3);
    if ~isempty(temp)
        if AMlbl_topology(temp(1),temp(end))==0
            AMlbl_topology(temp(1),temp(end))=L(i);
        else
            AMlbl_topology(temp(end),temp(1))=L(i);
        end
    end
end

% % cut and relabel AMlbl_topology
% cut=(sum(AMlbl_topology)==0);
% if nnz(cut)>0
%     AMlbl_topology(cut,:)=[];
%     AMlbl_topology(:,cut)=[];
%     AMlbl_topology = LabelBranchesAM(AMlbl_topology>0);
% end


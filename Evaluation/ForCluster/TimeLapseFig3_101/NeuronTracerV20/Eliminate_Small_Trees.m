% This function finds and eliminates small trees. 

function [AMlbl,r,R] = Eliminate_Small_Trees(AMlbl,r,R,L_thr)

disp('Eliminating small trees.')

TreeLabels=unique(AMlbl(AMlbl>0));
N_trees=length(TreeLabels);
for i=1:N_trees
    [ii,jj]=find(AMlbl==TreeLabels(i));
    L_tree=sum(sum((r(ii,:)-r(jj,:)).^2,2).^0.5)/2;
    if L_tree<=L_thr
        AMlbl(ii,jj)=0;
    end
end
cut=(sum(AMlbl)==0);
if nnz(cut)>0
    AMlbl(cut,:)=[];
    AMlbl(:,cut)=[];
    r(cut,:)=[];
    R(cut)=[];
    AMlbl = LabelTreesAM(AMlbl);
end

disp(['There are ',num2str(length(unique(AMlbl(AMlbl>0)))), ' trees in the stack.'])
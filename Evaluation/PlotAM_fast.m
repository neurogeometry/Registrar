% This version plots using fewer line objects (which slows down plotting when the number of vertices is large)
% The strategy is outlined below
% 1. Cycle through labeled trees
% 2. In each tree, first label unique branches
% 3. Order the vertices using AM
% 4. Plot each branch as a single line object

function []=PlotAM_fast(AM,r)
AM = max(AM,AM');
AM = triu(AM);

Labels=unique(full(AM(AM(:)>0)));
L=length(Labels);

if L==1
    colors=[1,0,0];
else
    colors=hsv(3*L);
    colors=colors(randperm(3*L),:);
    
    ind=(sum((colors-ones(size(colors,1),1)*[1,0,0]).^2,2).^0.5>0.3);
    colors=colors(ind,:);
    colors=colors(1:L,:);
    colors(1,:)=[1,0,0];
end

for f=1:L
    [i,j]=find(AM==Labels(f));
    ind=unique([i(:);j(:)]);
    AM_tree=AM(ind,:);AM_tree=AM_tree(:,ind);
    r_tree=r(ind,:);
    AMlbl_tree = LabelBranchesAM(AM_tree);
    BranchLabels=unique(AMlbl_tree(AMlbl_tree>0));
    for b=1:numel(BranchLabels)
        
        [i,j]=find(AMlbl_tree==BranchLabels(b));
        indb=unique([i(:);j(:)]);
        AM_branch=AMlbl_tree(indb,:);AM_branch=AM_branch(:,indb);
        r_branch=r_tree(indb,:);
        startt=find(sum(AM_branch>0,1)==1);
        [AM_branch,r_branch] = orderbranch(AM_branch,r_branch,startt(1));
        plot3(r_branch(:,2),r_branch(:,1),r_branch(:,3),'Color',colors(f,:),'LineStyle','-','LineWidth',1,'Visible','off');
    end
end
set(findobj(gca,'Type','Line'),'Visible','on')
end

function [AM,r] = orderbranch(AM,r,startt)
%Order r and AM. %AM labels are removed.
AM = max(AM,AM');
AM=spones(AM);
del=find(sum(AM,2)==0);
AM(del,:)=[];AM(:,del)=[];r(del,:)=[];

%Step along trace assuming no branching
ind=zeros(size(AM,1),1);
ind(1)=startt;
endpt=find(sum(AM,1)==1);
endpt(endpt==startt)=[];
[~,ind(2)]=find(AM(ind(1),:));
i=3;
while ind(i-1)~=endpt
    [~,next]=find(AM(ind(i-1),:));
    ind(i)=next(next~=ind(i-2));
    i=i+1;
end

AM=AM(ind,ind);
r=r(ind,:); 
end
% This function plots the tree structure contained in AM.
% The function works with labeled or not labeled AM.
% AM can be directed or undirected. 
% The labels don't have to be consecutive.

function h=PlotAM_h(AM,r)

AM = max(AM,AM');
AM = triu(AM);

Labels=unique(full(AM(AM(:)>0)));
L=length(Labels);

if L==1
    colors=[0.5,0.5,0.5]; %[1,0,0];
else
    colors=hsv(3*L);
    colors=colors(randperm(3*L),:);

    ind=(sum((colors-ones(size(colors,1),1)*[1,0,0]).^2,2).^0.5>0.3);
    colors=colors(ind,:);
    colors=colors(1:L,:);
    colors(1,:)=[1,0,0];
end

h=cell(1,L);
for f=1:L
    [i,j]=find(AM==Labels(f));
    X=[r(i,1),r(j,1)]';
    Y=[r(i,2),r(j,2)]';
    Z=[r(i,3),r(j,3)]';
    h{f}=line(Y,X,Z,'Color',colors(f,:),'LineStyle','-','LineWidth',3);
end
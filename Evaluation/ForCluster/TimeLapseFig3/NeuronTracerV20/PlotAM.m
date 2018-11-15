% This function plots the tree structure contained in AM.
% The function works with labeled or not labeled AM.
% AM can be directed or undirected. 
% The labels don't have to be consecutive.

function PlotAM(AM,r,color)

AM = max(AM,AM');
AM = triu(AM);

Labels=full(AM(AM(:)>0));
L_largest=mode(Labels);
Labels=unique(Labels);
L=length(Labels);

if L==1
    colors=rand(1,3); %[1,0,0]
else
    colors=hsv(3*L);
    colors=colors(randperm(3*L),:);
    
    ind=(sum((colors-ones(size(colors,1),1)*[1,0,0]).^2,2).^0.5>0.3);
    colors=colors(ind,:);
    colors=colors(1:L,:);
    colors(L_largest,:)=[1,0,0];
    %colors=parula(3*L);
    %colors=colors(randperm(size(colors,1)),:);
end

for f=1:L
    [i,j]=find(AM==Labels(f));
    X=[r(i,1),r(j,1)]';
    Y=[r(i,2),r(j,2)]';
    Z=[r(i,3),r(j,3)]';
    line(Y,X,Z,'Color',color,'LineWidth',1.5);
%     line(Y,X,Z,'Color',colors(f,:),'LineWidth',1);
end
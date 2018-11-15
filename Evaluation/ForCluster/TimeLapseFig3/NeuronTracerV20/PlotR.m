% This function plots the tree structures contained in AM.
% AM can be directed or undirected. 
% Brancg radii (R) are shown with colorbar

function PlotR(AM,r,R)

AM = spones(AM+AM');
AM = triu(AM);

[it,jt]=find(AM);
X=[r(it,1),r(jt,1)]';
Y=[r(it,2),r(jt,2)]';
Z=[r(it,3),r(jt,3)]';
R=(R(it)+R(jt))./2;

colors=jet(256);
C=colors(fix((R-min(R))./(max(R)-min(R)).*(256-1))+1,:);

for i=1:length(R)
    line(Y(:,i),X(:,i),Z(:,i),'Color',C(i,:),'LineWidth',2);
end
% This function plots the tree structure contained in AM.
% The function works with labeled or not labeled AM.
% AM can be directed or undirected. 
% The labels don't have to be consecutive.
% Color of plot is provided from outside.

function [h]=PlotAM_c(AM,r,color)

AM = max(AM,AM');
AM = triu(AM);

[i,j]=find(AM~=0);
X=[r(i,1),r(j,1)]';
Y=[r(i,2),r(j,2)]';
Z=[r(i,3),r(j,3)]';
h=line(Y,X,Z,'Color',color,'LineWidth',1);
h=h(1);

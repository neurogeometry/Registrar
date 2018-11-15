% This function plots the tree structure from an SWC format.
% Coordinates in the SWC file use Java frame of reference
% 0.5 is added to convert to Matlab format

function PlotSWC(pth_or_swc)

if ischar(pth_or_swc) && exist(pth_or_swc,'file')
    swc=textread(pth_or_swc);
elseif isnumeric(pth_or_swc) && size(pth_or_swc,2)==7
    swc=pth_or_swc;
else
    error('Incorrect path or SWC file format')
end

id=swc(:,1);
pid=swc(:,7);
r=[swc(:,4),swc(:,3),swc(:,5)]+0.5; % r ranges from 0.5 to sizeIm+0.5 in Matlab and [0 sizeIm] in Java and SWC

Roots=find(pid==-1);
plot3(r(Roots,2),r(Roots,1),r(Roots,3),'m.','MarkerSize',10)
hold on

for i=[id(pid~=-1)]'
    line([r(pid(i),2),r(i,2)],[r(pid(i),1),r(i,1)],[r(pid(i),3),r(i,3)],'Color','m')
end

% This function plots the tree structures (in different colors) from an SWC format.
% Coordinates in the SWC file use Java frame of reference
% 0.5 is added to convert to Matlab format

function PlotSWCtrees(pth_or_swc)

if exist(pth_or_swc,'file')
    swc=textread(pth_or_swc);
elseif size(pth_or_swc,2)==7
    swc=pth_or_swc;
else
    error('Incorrect path or SWC file format')
end

id=swc(:,1);
pid=swc(:,7);
r=[swc(:,4),swc(:,3),swc(:,5)]+0.5; % r ranges from 0.5 to sizeIm+0.5 in Matlab and [0 sizeIm] in Java and SWC

Roots=find(pid==-1);

if length(Roots)==1
    colors='g';
else
    colors=hsv(length(Roots));
    colors=colors(randperm(length(Roots)),:);
end

hold on
for i=1:length(Roots)
    kids=zeros(1,length(id));
    currentV=find(pid==Roots(i));
    while ~isempty(currentV)
        kids(nnz(kids)+1:nnz(kids)+length(currentV))=currentV;
        [currentV,temp]=find(ones(length(pid),1)*currentV'-pid*ones(1,length(currentV))==0);
    end
    kids(kids==0)=[];
    plot3(r(Roots(i),2),r(Roots(i),1),r(Roots(i),3),'Marker','.','Color',colors(i,:),'MarkerSize',10)
    line([r(pid(kids),2),r(kids,2)]',[r(pid(kids),1),r(kids,1)]',[r(pid(kids),3),r(kids,3)]','Color',colors(i,:))
end



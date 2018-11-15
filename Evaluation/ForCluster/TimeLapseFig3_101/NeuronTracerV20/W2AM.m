% This function converts the wave trace (sparse format) in an AM r format
% Works also for sparse Im 

function [AM,r]=W2AM(Im,W,sizeIm)

if isempty(sizeIm)
    Im=double(Im);
    sizeIm=size(Im);
    if length(sizeIm)==2
        sizeIm(3)=1;
    end
end

sizeAM=max(W);
AM=sparse(sizeAM,sizeAM); 
r=zeros(sizeAM,3);

SE=strel(ones(3,3,3));
offsets=getneighbors(SE);
ii=offsets(:,1); jj=offsets(:,2); kk=offsets(:,3);

[ind,~,Labels]=find(W);
[x,y,z]=ind2sub_AS(sizeIm,find(W(:)==1,1,'first'));
r(1,:)=[x y z];

for i=1:sizeAM 
    NextLabels=unique(W(RegionNhood(ind(find(Labels==i)),sizeIm,ii,jj,kk)));
    NextLabels(NextLabels==0)=[];
    NextLabels(NextLabels==i)=[];
    
    for j=1:length(NextLabels)
        Label_ind=ind(find(Labels==NextLabels(j))); % find is faster here
        [x,y,z]=ind2sub_AS(sizeIm,Label_ind);
        r(NextLabels(j),:)=Im(Label_ind)'*[x,y,z]./(sum(Im(Label_ind)));
        AM(i,NextLabels(j))=1;
    end
end
AM=logical(AM+AM');
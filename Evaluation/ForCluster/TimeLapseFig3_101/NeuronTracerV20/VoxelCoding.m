% This function propagates a wave of labels through the connected regions in
% the image. The wave originates from the start voxel (if provided). 
% If flag==0 the wave propagates into only one region connected to the start voxel.
% W is returned in sparse format

function W=VoxelCoding(Im,Start_X,Start_Y,Start_Z,flag,varargin)

if length(varargin)==1
    myClass=varargin{1};
else
    myClass.getFlag=NaN;
end

sizeIm=size(Im);
if length(sizeIm)==2
    sizeIm(3)=1;
end

StartVoxel=sub2ind_AS(sizeIm,Start_X,Start_Y,Start_Z);
if isempty(StartVoxel)
    StartVoxel=find(Im,1);
end

Front=StartVoxel;

SE=strel(ones(3,3,3));
[offsets, ~]=getneighbors(SE);
ii=offsets(:,1); jj=offsets(:,2); kk=offsets(:,3);

W=sparse(prod(sizeIm),1);
W(Front)=1; 
AvlFronts={}; 
n=1;
NewFront=1;
stop_cond=1;

while stop_cond
    if myClass.getFlag()==1
        return;
    end
    
    NHood=RegionNhood(Front,sizeIm,ii,jj,kk);
    NHood=NHood(W(NHood)==0 & Im(NHood)~=0);
        
    if isempty(NHood)
        if ~isempty(AvlFronts)
            NewFront=AvlFronts{end};
            NewFront=NewFront(W(NewFront)==0);
            AvlFronts(end)=[];
        else
            NewFront=[];
        end
        
    else
        [NbrFronts NF]=ChkBranching_AS(NHood,sizeIm);
        if NbrFronts==1
            NewFront=NF{1}; 
        elseif NbrFronts>1
            for f=1:length(NF)-1
                ix = length(AvlFronts)+1;
                AvlFronts(ix)=NF(f);
            end
            NewFront=NF{end};  
        end
    end
    
    if ~isempty(NewFront)
        n = n + 1;
        W(NewFront)=n;
        Front=NewFront;
    end
    
    stop_cond=(length(AvlFronts)>1 || ~isempty(NewFront));
    if stop_cond==0 && (isempty(Start_X+Start_Y+Start_Z) || flag==1)
        Front=find((W>0)~=(Im(:)>0),1);
        if ~isempty(Front)
            stop_cond=1;
            n=n+1;
            W(Front)=n;
        end
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function determines the indexes of all connected regions in the
% branched front. Cell array NF contains these indexes
function [NbrFronts NF] = ChkBranching_AS(NewFront,sizeIm)

NbrFronts=0;
NF=[];

if ~isempty(NewFront)
    [x y z]=ind2sub_AS(sizeIm,NewFront);
    MIx=min(x); MIy=min(y); MIz=min(z);
    MAx=max(x); MAy=max(y); MAz=max(z);
    
    sizeIm_small=[MAx-MIx+1,MAy-MIy+1,MAz-MIz+1];
    ind_small=sub2ind_AS(sizeIm_small,x-MIx+1,y-MIy+1,z-MIz+1);
    Im_small=zeros(sizeIm_small);
    Im_small(ind_small)=1;
    
    zx=bwlabeln(Im_small);
    NbrFronts=max(zx(:));
    
    NF=cell(1,NbrFronts);
    if NbrFronts>1
        for f=1:NbrFronts
            rr=find(zx(:)==f);
            [x1 y1 z1]=ind2sub_AS(sizeIm_small,rr);
            NF{f}=sub2ind_AS(sizeIm,MIx+x1-1,MIy+y1-1,MIz+z1-1);
        end
    else
        NF{1}=NewFront;
    end
end  
    
    
    
    
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% this function generates N uniformly distributed seed points
% Mesh is the minimum distance between the seeds
% thr is the threshold
% R_min,R_step, and R_max are parameters of Multi_Scale_LoG filter

function SVr=Find_Seeds(Orig,Mesh,thr,R_min,R_step,R_max,SVr,varargin)

if length(varargin)==1
    myClass=varargin{1};
else
    myClass.getFlag=NaN;
end

Orig=Multi_Scale_LoG(Orig,R_min,R_step,R_max,myClass);

sizeIm=size(Orig);
if length(sizeIm)==2
    sizeIm=[sizeIm,1];
end

for i=1:size(SVr,1)
    if myClass.getFlag()==1
        return;
    end
    Orig(max(SVr(i,1)-Mesh,1):min(SVr(i,1)+Mesh,sizeIm(1)),max(SVr(i,2)-Mesh,1):min(SVr(i,2)+Mesh,sizeIm(2)),max(SVr(i,3)-Mesh,1):min(SVr(i,3)+Mesh,sizeIm(3)))=0;
end

N=ceil(sizeIm(1)/Mesh)*ceil(sizeIm(2)/Mesh)*ceil(sizeIm(3)/Mesh);
SVx=zeros(N,1);
SVy=zeros(N,1);
SVz=zeros(N,1);

count=0;
[Max,ind]=max(Orig(:));
while Max>=thr
    if myClass.getFlag()==1
        return;
    end
    count=count+1;
    [SVx(count),SVy(count),SVz(count)]=ind2sub(sizeIm,ind);
    Orig(max(SVx(count)-Mesh,1):min(SVx(count)+Mesh,sizeIm(1)),max(SVy(count)-Mesh,1):min(SVy(count)+Mesh,sizeIm(2)),max(SVz(count)-Mesh,1):min(SVz(count)+Mesh,sizeIm(3)))=0;
    [Max,ind]=max(Orig(:));
end

SVr=[SVr;[SVx(1:count),SVy(1:count),SVz(1:count)]];

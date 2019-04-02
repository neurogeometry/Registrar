%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% this function generates N uniformly distributed seed points
% Mesh is the minimum distance between the seeds
% thr is the threshold
% R_min,R_step, and R_max are parameters of Multi_Scale_LoG filter

function SVr=Find_Seeds(Orig,FEsmooth3BoxSize,FEk,FEExpected_Missalignment,FEMaxN_Features,FEMesh)
% FEsmooth3BoxSize = 3;
% FEk = 0.5;
% FEExpected_Missalignment = [0.1,0.1,0.1];
% FEMaxN_Features=1000;
% FEMesh = [15 15 5];

Orig = double(Orig);
SizeIM = size(Orig);

Orig=smooth3(Orig,'box',FEsmooth3BoxSize);

thr = mean(Orig(:))+FEk*std(Orig(:));

sizeIm=size(Orig);
if length(sizeIm)==2
    sizeIm=[sizeIm,1];
end

ind=find(Orig>thr);
keep_ind=false(1,length(ind));
[x,y,z]=ind2sub(sizeIm,ind);
delt=round(FEExpected_Missalignment.*sizeIm);

% if length (SizeIM) ~= 2
%     for i=2:size(StackPositions_pixels,1)
%         keep_ind(x+StackPositions_pixels(1,1)-1+delt(1)>=StackPositions_pixels(i,1) & x+StackPositions_pixels(1,1)-1-delt(1)<=StackPositions_pixels(i,1)+StackSizes_pixels(i,1)-1 &...
%             y+StackPositions_pixels(1,2)-1+delt(2)>=StackPositions_pixels(i,2) & y+StackPositions_pixels(1,2)-1-delt(2)<=StackPositions_pixels(i,2)+StackSizes_pixels(i,2)-1 &...
%             z+StackPositions_pixels(1,3)-1+delt(3)>=StackPositions_pixels(i,3) & z+StackPositions_pixels(1,3)-1-delt(3)<=StackPositions_pixels(i,3)+StackSizes_pixels(i,3)-1)=true;
%     end
%     ind=ind(keep_ind);
%     x=x(keep_ind);
%     y=y(keep_ind);
%     z=z(keep_ind);
% end

clear keep_ind
Orig=Orig(ind);

[~,ind]=sort(Orig,'descend');
clear Orig
x=x(ind);
y=y(ind);
z=z(ind);

count=0;
temp_ind=1;
SVx=zeros(2*FEMaxN_Features,1);
SVy=zeros(2*FEMaxN_Features,1);
SVz=zeros(2*FEMaxN_Features,1);

while ~isempty(temp_ind) && count<=FEMaxN_Features
    count=count+1;
%     disp(x);
    SVx(count)=x(temp_ind);
    SVy(count)=y(temp_ind);
    SVz(count)=z(temp_ind);
    rem_ind=find(abs(x-x(temp_ind))<=FEMesh(1) & abs(y-y(temp_ind))<=FEMesh(2) & abs(z-z(temp_ind))<=FEMesh(3));
    ind(rem_ind)=0;
    temp_ind=find(ind>0,1,'first');
end

SVr=[SVx(1:count),SVy(1:count),SVz(1:count)];

% This function solves of the Eikonal equation by using the Fast Marching algorithm 
% of Sethian. T is the arival time map. The boundary condition is set as 
% T(Start_X,Start_Y,Start_Z) = 0 (if start point is provided). 

function T=FastMarching_old(Orig)

disp('Fast Marching started.')

N_steps=200; % number of re-initialization steps
Max_dist=15; % re-initialization is done if a fornt reaches Max_dist

Orig=Multi_Scale_LoG(Orig,2,1,3);
Orig=double(Orig);
Orig=Orig-min(Orig(:));
Orig=(Orig./max(Orig(:)));
Orig(Orig<0.01)=0.01;

sizeOrig=size(Orig);
if length(sizeOrig)==2
    sizeOrig=[sizeOrig,1];
end
Im=zeros(sizeOrig+4);
Im(3:end-2,3:end-2,3:end-2)=Orig;
sizeIm=size(Im);

BW = imregionalmax(Im,26);
StartVoxel=find(BW & Im>0.75);
StartVoxel=StartVoxel(10);
[SVx,SVy,SVz]=ind2sub(sizeIm,StartVoxel);

SE=zeros(3,3,3);
SE([5,11,13,15,17,23])=1; % 6-neighbors
SE=strel(SE);
[offsets, ~]=getneighbors(SE);
ii=offsets(:,1); jj=offsets(:,2); kk=offsets(:,3);
iip=ii'; jjp=jj'; kkp=kk';

T=inf(sizeIm); T(StartVoxel)=0;
D=inf(sizeIm); D(StartVoxel)=0;
L=zeros(sizeIm); L(StartVoxel)=(1:length(StartVoxel));

KnownPoints=zeros(10000,1);
NPcount=length(StartVoxel);
KnownPoints(1:NPcount)=StartVoxel;
TrialPoints=zeros(10000,1);
TPcount=0;
T_TrialPoints=inf(10000,1);

NHood=RegionNhood(KnownPoints(1:NPcount),sizeIm,ii,jj,kk);
NHood(ismember(NHood,KnownPoints(1:NPcount)))=[];

figure(100)
imshow(max(Im,[],3)), hold on
[NVx,NVy,NVz]=ind2sub(sizeIm,KnownPoints(1:NPcount));
plot3(NVy,NVx,NVz,'r.')
drawnow

for i=1:N_steps 
    display([i/N_steps,length(unique(L(L>0)))]);
    stop_cond=true;
    while stop_cond      
%         [x,y,z]=ind2sub_AS(sizeIm,NHood);
%         X=ones(length(x),1)*iip+x*ones(1,6);
%         Y=ones(length(y),1)*jjp+y*ones(1,6);
%         Z=ones(length(z),1)*kkp+z*ones(1,6);
%         ind = X+(Y-1).*sizeIm(1)+(Z-1).*(sizeIm(1)*sizeIm(2));
        ind=[NHood-1,NHood+1,NHood-sizeIm(1),NHood+sizeIm(1),NHood-sizeIm(1)*sizeIm(2),NHood+sizeIm(1)*sizeIm(2)];
        [T(NHood),D(NHood)]=BackSolver(Im(NHood),T(ind),D(ind)); 
        keep_ind=(L(NHood)==0);
        L(NHood(keep_ind))=max(L(ind(keep_ind,:)),[],2);
        TrialPoints(TPcount+1:TPcount+sum(keep_ind))=NHood(keep_ind);
        T_TrialPoints(TPcount+1:TPcount+sum(keep_ind))=T(NHood(keep_ind));
        TPcount=TPcount+sum(keep_ind);
        
        [~,ind2]=min(T_TrialPoints); %(1:TPcount)
        NewKnownPoint=TrialPoints(ind2);
        TrialPoints(ind2)=inf;
        T_TrialPoints(ind2)=inf;
        KnownPoints(NPcount+1)=NewKnownPoint;
        NPcount=NPcount+1;
        [x,y,z]=ind2sub_AS(sizeIm,NewKnownPoint);
        NHood = ii+x+(jj+y-1).*sizeIm(1)+(kk+z-1).*(sizeIm(1)*sizeIm(2));
        collision=(sum(L(NHood))~=sum(L(NHood)>0).*max(L(NHood)));
        %NHood(ismember(NHood,KnownPoints))=[];
        NHood(T(NHood)<=T(NewKnownPoint))=[];
                
        if D(NewKnownPoint)>=Max_dist
            stop_cond=false;
        end
        if collision
            stop_cond=false;
        end
    end
    
    [T,D,L]=GradientDescent(T,D,L,NewKnownPoint);
    ind=find(TrialPoints(1:TPcount)==inf); % remove inf
    TrialPoints(ind)=[];
    T_TrialPoints(ind)=[];
    TPcount=TPcount-length(ind);
    
    ind=find(T(TrialPoints(1:TPcount))==inf | L(TrialPoints(1:TPcount))==L(NewKnownPoint));
    TrialPoints(ind)=[];
    T_TrialPoints(ind)=[];
    TPcount=TPcount-length(ind);
    
    ind=find(T(KnownPoints(1:NPcount))==inf | L(KnownPoints(1:NPcount))==L(NewKnownPoint));
    KnownPoints(ind)=[];
    NPcount=NPcount-length(ind);
    new_region_ind=find(L==L(NewKnownPoint));
    KnownPoints(NPcount+1:NPcount+length(new_region_ind))=new_region_ind;
    NPcount=NPcount+length(new_region_ind);
    
    NHood=RegionNhood(new_region_ind,sizeIm,ii,jj,kk);
    NHood(ismember(NHood,KnownPoints(1:NPcount)))=[];
    
    figure(100)
    [NVx,NVy,NVz]=ind2sub(sizeIm,new_region_ind);
    plot3(NVy,NVx,NVz,'r.')
    drawnow
end

Labels=unique(L(:));
Labels(Labels==0)=[];
for i=1:length(Labels)
    temp_ind=(L==Labels(i));
    temp_ind2=(T>min(T(temp_ind)));
    
    L(temp_ind & temp_ind2)=0;
    D(temp_ind & temp_ind2)=inf;
    T(temp_ind & temp_ind2)=inf;
end

disp('Fast Marching is complete.')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% [EVx,EVy,EVz]=ind2sub(sizeIm,NewKnownPoint);
% 
figure
imshow(max(L,[],3),[0 max(L(:))]), hold on
plot3(SVy,SVx,SVz,'r.')
plot3(EVy,EVx,EVz,'g.')
% 
% figure
% Ttemp=T;
% Ttemp(T==inf)=0;
% imshow(max(Ttemp,[],3),[0 max(Ttemp(:))]), hold on
% plot3(SVy,SVx,SVz,'r.')
% plot3(EVy,EVx,EVz,'g.')
% 
% figure
% Dtemp=D;
% Dtemp(D==inf)=0;
% Dtemp(D==0)=max(Dtemp(:));
% imshow(max(Dtemp,[],3),[0 max(Dtemp(:))]), hold on
% plot3(SVy,SVx,SVz,'r.')
% plot3(EVy,EVx,EVz,'g.')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% this function finds arrival time, T, and distance, D, for the trial points
function [T,D]=BackSolver(F,Tpm_xyz,Dpm_xyz)

delta=1;
N=ones(size(Tpm_xyz,1),1)*(1:size(Tpm_xyz,2));
ct=(delta./F).^2;
cd=(delta).^2.*ones(size(F));
Ct=ct*ones(1,size(Tpm_xyz,2));
Cd=cd*ones(1,size(Dpm_xyz,2));
Tpm_xyz=sort(Tpm_xyz,2);
Dpm_xyz=sort(Dpm_xyz,2);
indt=((cumsum(Tpm_xyz.^2,2)-cumsum(Tpm_xyz,2).^2./N)<=Ct);
indd=((cumsum(Dpm_xyz.^2,2)-cumsum(Dpm_xyz,2).^2./N)<Cd); % CHECK
Tpm_xyz(~indt)=0;
Dpm_xyz(~indd)=0;
nt=sum(indt,2);
nd=sum(indd,2);
Tmean=sum(Tpm_xyz,2)./nt;
Dmean=sum(Dpm_xyz,2)./nd;
T2mean=sum(Tpm_xyz.^2,2)./nt;
D2mean=sum(Dpm_xyz.^2,2)./nd;
T=Tmean+(ct./nt-(T2mean-Tmean.^2)).^0.5;
D=Dmean+(cd./nd-(D2mean-Dmean.^2)).^0.5;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function performs gradient descent from ind (a single point) to
% T=0. T map of the region is reset to the gradient descent path, Tnew.
% If there are multiple labels neighboring ind, the gradient descent will
% go into all these labeld regions. The result will be merged and
% relabeled.
function [T,D,L]=GradientDescent(T,D,L,ind)

sizeIm=size(T);
T_ind=T(ind);

SE=zeros(3,3,3);
SE([5,11,13,15,17,23])=1; % 6-neighbors
SE=strel(SE);
[offsets, ~]=getneighbors(SE);
ii=offsets(:,1); jj=offsets(:,2); kk=offsets(:,3);

%[x,y,z]=ind2sub_AS(sizeIm,ind);
%pot_ind = (ii+x)+(jj+y-1).*sizeIm(1)+(kk+z-1).*(sizeIm(1)*sizeIm(2));
pot_ind=[ind-1,ind+1,ind-sizeIm(1),ind+sizeIm(1),ind-sizeIm(1)*sizeIm(2),ind+sizeIm(1)*sizeIm(2)];
Labels=unique(L(pot_ind));
Labels(Labels==0)=[];

if length(Labels)>1
    L_new=max(L(:))+1;
else
    L_new=L(ind);
end

T_0=zeros(1,length(Labels));
for i=1:length(Labels)
    T_0(i)=min(T(L==Labels(i)));
end

all_ind=zeros(1,sum(sizeIm));
all_ind(1)=ind;
count=1;
for i=1:length(Labels)
    current_ind=ind;
    while T(current_ind)~=T_0(i) && count<1000 %%% FIX
        %[x,y,z]=ind2sub_AS(sizeIm,current_ind);
        %pot_ind = (ii+x)+(jj+y-1).*sizeIm(1)+(kk+z-1).*(sizeIm(1)*sizeIm(2));
        pot_ind=[current_ind-1,current_ind+1,current_ind-sizeIm(1),current_ind+sizeIm(1),current_ind-sizeIm(1)*sizeIm(2),current_ind+sizeIm(1)*sizeIm(2)];
        pot_ind(L(pot_ind)~=Labels(i))=[];
        
        [~,temp]=min(T(pot_ind));
        current_ind=pot_ind(temp);
        count=count+1;
        if count==1000
            disp(count)
        end
        all_ind(count)=current_ind;
    end
end
all_ind=all_ind(1:count);

for i=1:length(Labels)
    temp_ind=(L==Labels(i));
    temp_ind1=(T>T_0(i));
    temp_ind2=(T==T_0(i));
    
    L(temp_ind & temp_ind1)=0;
    L(temp_ind & temp_ind2)=L_new;
    
    D(temp_ind & temp_ind1)=inf;
    D(temp_ind & temp_ind2)=0;
    
    T(temp_ind & temp_ind1)=inf;
    T(temp_ind & temp_ind2)=T_ind;
end

L(all_ind)=L_new;
D(all_ind)=0;
T(all_ind)=T_ind;




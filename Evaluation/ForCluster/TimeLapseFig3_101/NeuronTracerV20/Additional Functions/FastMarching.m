% This function solves of the Eikonal equation by using the Fast Marching algorithm 
% of Sethian. T is the arival time map. The boundary condition is set as 
% T(Start_X,Start_Y,Start_Z) = 0 (if start point is provided). 

function T=FastMarching(Orig)

disp('Fast Marching started.')

N_steps=100; % number of re-initialization steps
Max_dist=10; % re-initialization is done if a fornt reaches Max_dist
seed_R=5; 
seed_thr=80;

Orig=Multi_Scale_LoG(Orig,2,1,2);

[SVx,SVy,SVz,~]=Find_Seeds(Orig,seed_R,seed_thr);

sizeOrig=size(Orig);
if length(sizeOrig)==2
    sizeOrig=[sizeOrig,1];
end

Orig=double(Orig);
%Orig=Orig-min(Orig(:));
%Orig=(Orig./max(Orig(:))).^1;
%Orig(Orig<0.01)=0.01;

pad=2;
Im=zeros(sizeOrig+2*pad);
Im(1+pad:end-pad,1+pad:end-pad,1+pad:end-pad)=Orig;
sizeIm=size(Im);
SVx=SVx+pad;
SVy=SVy+pad;
SVz=SVz+pad;
StartVoxel=sub2ind_AS(sizeIm,SVx,SVy,SVz);
StartVoxel=StartVoxel(2);

N6_ind=[-1;+1;-sizeIm(1);+sizeIm(1);-sizeIm(1)*sizeIm(2);+sizeIm(1)*sizeIm(2)];

T=inf(sizeIm); T(StartVoxel)=0;
D=inf(sizeIm); D(StartVoxel)=0;
L=zeros(sizeIm); L(StartVoxel)=(1:length(StartVoxel));
KT=zeros(sizeIm); KT(StartVoxel)=1; % Known=1, Trial=2

KnownPoints=zeros(10000,1);
NPcount=length(StartVoxel);
KnownPoints(1:NPcount)=StartVoxel;
TrialPoints=zeros(10000,1);
TPcount=0;
T_TrialPoints=inf(10000,1);

NHood=[KnownPoints(1:NPcount)+N6_ind(1);KnownPoints(1:NPcount)+N6_ind(2);KnownPoints(1:NPcount)+N6_ind(3); ...
    KnownPoints(1:NPcount)+N6_ind(4);KnownPoints(1:NPcount)+N6_ind(5);KnownPoints(1:NPcount)+N6_ind(6)];
NHood=unique(NHood);
NHood(KT(NHood)==1)=[];

figure(100)
imshow(max(Im,[],3),[0 max(Im(:))]), hold on
[NVx,NVy,NVz]=ind2sub(sizeIm,KnownPoints(1:NPcount));
plot3(NVy,NVx,NVz,'g.')
drawnow

%tic
for i=1:N_steps 
    display([i/N_steps,length(unique(L(KnownPoints(1:NPcount))))]);
    stop_cond=true;
    while stop_cond      
        ind=[NHood+N6_ind(1),NHood+N6_ind(2),NHood+N6_ind(3),NHood+N6_ind(4),NHood+N6_ind(5),NHood+N6_ind(6)];
       
        %[T(NHood),D(NHood)]=BackSolver(Im(NHood),T(ind),D(ind));
        %[~,temp_ind]=min(T(ind),[],2);
        
        temp_T=T(ind);
        temp_T(KT(ind)~=1)=inf;
        [T(NHood),D(NHood)]=BackSolver(Im(NHood),temp_T,D(ind));
        
        [~,temp_ind]=min(temp_T,[],2);
        
        temp=(1:length(NHood))'+length(NHood).*(temp_ind-1);
        L(NHood)=L(ind(temp));
        
        keep_ind=(KT(NHood)==0);
        %KT(NHood)=2;
        KT(NHood(keep_ind))=2;
        TrialPoints(TPcount+1:TPcount+sum(keep_ind))=NHood(keep_ind);
        T_TrialPoints(TPcount+1:TPcount+sum(keep_ind))=T(NHood(keep_ind));
        TPcount=TPcount+sum(keep_ind);
        
        [~,ind2]=min(T_TrialPoints); %(1:TPcount)
        NewKnownPoint=TrialPoints(ind2);
        KT(NewKnownPoint)=1;
        TrialPoints(ind2)=inf;
        T_TrialPoints(ind2)=inf;
        KnownPoints(NPcount+1)=NewKnownPoint;
        NPcount=NPcount+1;
        NHood=NewKnownPoint+N6_ind;
        
        tempL=L(NHood(KT(NHood)==1)); %L(NHood);
        tempL=tempL(tempL>0);
        collision=(sum(tempL-tempL(1))~=0);
        NHood(KT(NHood)==1)=[];
        %KT(NHood)=2;

        if D(NewKnownPoint)>=Max_dist
            stop_cond=false;
        end
        if collision
            stop_cond=false;
        end
    end
    
    [T,D,L,KT]=GradientDescent(T,D,L,KT,NewKnownPoint);
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
    
    NHood=[new_region_ind+N6_ind(1);new_region_ind+N6_ind(2);new_region_ind+N6_ind(3); ...
        new_region_ind+N6_ind(4);new_region_ind+N6_ind(5);new_region_ind+N6_ind(6)];
    NHood=unique(NHood);
    NHood(KT(NHood)==1)=[];
    
%     figure(200),hold on
%     plot(NPcount,toc,'g*')
%     drawnow
    
    figure(100)
    [NVx,NVy,NVz]=ind2sub(sizeIm,new_region_ind);
    plot3(NVy,NVx,NVz,'r.')
    drawnow
end

Labels=unique(L(:));
Labels(Labels==0)=[];
for i=1:length(Labels)
    temp_ind1=(L==Labels(i));
    temp_ind2=(T>min(T(temp_ind1)));
    
    L(temp_ind1 & temp_ind2)=0;
    D(temp_ind1 & temp_ind2)=inf;
    T(temp_ind1 & temp_ind2)=inf;
end

disp('Fast Marching is complete.')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[EVx,EVy,EVz]=ind2sub(sizeIm,NewKnownPoint);

figure
imshow(max(L,[],3),[0 max(L(:))]), hold on
plot3(SVy,SVx,SVz,'r.')
plot3(EVy,EVx,EVz,'g.')

figure
Ltemp=L;
Ltemp(KT==2)=0;
imshow(max(Ltemp,[],3),[0 max(Ltemp(:))]), hold on
plot3(SVy,SVx,SVz,'r.')
plot3(EVy,EVx,EVz,'g.')

figure
Ttemp=T;
Ttemp(T==inf)=0;
imshow(max(Ttemp,[],3),[0 max(Ttemp(:))]), hold on
plot3(SVy,SVx,SVz,'r.')
plot3(EVy,EVx,EVz,'g.')

figure
Dtemp=D;
Dtemp(D==inf)=0;
Dtemp(D==0)=max(Dtemp(:));
imshow(max(Dtemp,[],3),[0 max(Dtemp(:))]), hold on
plot3(SVy,SVx,SVz,'r.')
plot3(EVy,EVx,EVz,'g.')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% this function finds N evenly distributed seed points
function [SVx,SVy,SVz,StartVoxel]=Find_Seeds(Im,R,thr)

sizeIm=size(Im);
if length(sizeIm)==2
   sizeIm=[sizeIm,1]; 
end

N=10000;
StartVoxel=zeros(N,1);
SVx=zeros(N,1);
SVy=zeros(N,1);
SVz=zeros(N,1);
count=0; Max=thr;
while Max>=thr
    count=count+1;
    [Max,ind]=max(Im(:));
    StartVoxel(count)=ind;
    [SVx(count),SVy(count),SVz(count)]=ind2sub(sizeIm,StartVoxel(count));
    Im(max(SVx(count)-R,1):min(SVx(count)+R,sizeIm(1)),max(SVy(count)-R,1):min(SVy(count)+R,sizeIm(2)),max(SVz(count)-R,1):min(SVz(count)+R,sizeIm(3)))=0;
end

StartVoxel=StartVoxel(1:count-1);
SVx=SVx(1:count-1);
SVy=SVy(1:count-1);
SVz=SVz(1:count-1);

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
function [T,D,L,KT]=GradientDescent(T,D,L,KT,ind)

% if ind==778782
%     1
% end
sizeIm=size(T);
% T_ind=T(ind);
% reduce=0.5;
% D_reset=1;

N6_ind=[-1,+1,-sizeIm(1),+sizeIm(1),-sizeIm(1)*sizeIm(2),+sizeIm(1)*sizeIm(2)];
pot_ind=ind+N6_ind;
pot_ind(KT(pot_ind)~=1)=[];
Labels=unique(L(pot_ind));

if length(Labels)>1
    L_new=max(L(:))+1;
else
    L_new=L(ind);
end

T_0=zeros(1,length(Labels));
% for i=1:length(Labels)
%     T_0(i)=min(T(L==Labels(i)));
% end

all_ind=zeros(1,sum(sizeIm));
all_ind(1)=ind;
count=1;
for i=1:length(Labels)
    current_ind=ind;
    while T(current_ind)~=T_0(i) %&& count<100 %%% FIX
        pot_ind=current_ind+N6_ind;
        pot_ind(L(pot_ind)~=Labels(i))=[];
        pot_ind(KT(pot_ind)~=1)=[];
        pot_ind(pot_ind==ind)=[];
        
        [~,temp]=min(T(pot_ind));
        current_ind=pot_ind(temp);
        if ~isempty(current_ind)
            count=count+1;
            if count==100
                error(num2str(count))
            end
            all_ind(count)=current_ind;
        end
    end
end
all_ind=all_ind(1:count);

for i=1:length(Labels)
    temp_ind=find(L==Labels(i));
    temp_ind1=temp_ind(T(temp_ind)>T_0(i));
    temp_ind2=temp_ind(T(temp_ind)<=T_0(i));
    %temp_ind3=temp_ind(KT(temp_ind)==2);
    
%     temp_ind1=temp_ind(T(temp_ind)>T_ind.*reduce);
%     temp_ind2=temp_ind(T(temp_ind)<=T_ind.*reduce);
%     temp_ind1=temp_ind(D(temp_ind)>D_reset | T(temp_ind)==inf);
%     temp_ind2=temp_ind(D(temp_ind)<=D_reset & T(temp_ind)<inf);
    
    L(temp_ind1)=0;
    L(temp_ind2)=L_new;
    %L(temp_ind3)=0;
    
    D(temp_ind1)=inf;
    D(temp_ind2)=0;
    %D(temp_ind3)=inf;
    
    T(temp_ind1)=inf;
    T(temp_ind2)=0; %T_ind;
    %T(temp_ind3)=inf;
    
    KT(temp_ind1)=0;
    KT(temp_ind2)=1;
    %KT(temp_ind3)=0;
end

L(all_ind)=L_new;
D(all_ind)=0;
T(all_ind)=0; %T_ind;
KT(all_ind)=1;



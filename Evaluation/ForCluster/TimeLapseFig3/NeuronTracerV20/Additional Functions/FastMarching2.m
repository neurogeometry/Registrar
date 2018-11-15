% This function solves of the Eikonal equation by using the Fast Marching algorithm 
% of Sethian. T is the arival time map. The boundary condition is set as 
% T(Start_X,Start_Y,Start_Z) = 0 (if start point is provided). 

function L=FastMarching2(Orig)

disp('Fast Marching started.')

N_steps=900; % number of re-initialization steps
Max_dist=10; % re-initialization is done if a fornt reaches Max_dist
Max_dist=ceil(Max_dist);
seed_R=25; 
seed_thr=80;

Orig=Multi_Scale_LoG(Orig,2,1,4);

[SVx,SVy,SVz,~]=Find_Seeds(Orig,seed_R,seed_thr);

sizeOrig=size(Orig);
if length(sizeOrig)==2
    sizeOrig=[sizeOrig,1];
end

Orig=double(Orig);
Orig=Orig-min(Orig(:));
Orig=(Orig./max(Orig(:)));
Orig(Orig<0.01)=0.01;

pad=2;
Im=zeros(sizeOrig+2*pad);
Im(1+pad:end-pad,1+pad:end-pad,1+pad:end-pad)=Orig;
sizeIm=size(Im);
SVx=SVx+pad;
SVy=SVy+pad;
SVz=SVz+pad;
StartVoxel=sub2ind_AS(sizeIm,SVx,SVy,SVz);
%StartVoxel=StartVoxel(2);

N6_ind=[-1;+1;-sizeIm(1);+sizeIm(1);-sizeIm(1)*sizeIm(2);+sizeIm(1)*sizeIm(2)];

pad1=1*Max_dist+pad;
%sizeIm1=Max_dist+2*pad1.*[1,1,1];
%N6_ind1=[-1;+1;-sizeIm1(1);+sizeIm1(1);-sizeIm1(1)*sizeIm1(2);+sizeIm1(1)*sizeIm1(2)];

T=inf(sizeIm); T(StartVoxel)=0;
D=inf(sizeIm); D(StartVoxel)=0;
L=zeros(sizeIm); L(StartVoxel)=(1:length(StartVoxel));
KT=zeros(sizeIm); KT(StartVoxel)=1; % Known=1, Trial=2

figure(100)
imshow(max(Im,[],3),[0 max(Im(:))]), hold on
[NVx,NVy,NVz]=ind2sub(sizeIm,StartVoxel);
plot3(NVy,NVx,NVz,'g.')
drawnow

stop_condition=false;
count=0;
while stop_condition==false && count<=N_steps
    count=count+1;
    display([count,length(unique(L(KT==1)))]);

    % do fast marching untill 1 of 4 conditions is encountered
    [T,D,L,KT,NewKnownPoint,exit_flag]=Fast_Marching(Im,T,D,L,KT,N6_ind,Max_dist,inf);
    Max_T=T(NewKnownPoint);
    if strcmp(exit_flag,'dist') || strcmp(exit_flag,'collision')
        % find the gradient descent path and relabel the regions
        [T,D,L,KT,GDpath]=GradientDescent(T,D,L,KT,NewKnownPoint,N6_ind);
                
        [GDx,GDy,GDz]=ind2sub(sizeIm,GDpath);
        min_x=min(GDx); max_x=max(GDx);
        min_y=min(GDy); max_y=max(GDy);
        min_z=min(GDz); max_z=max(GDz);
        
        sizeIm1=[max_x-min_x,max_y-min_y,max_z-min_z]+1+2*pad1;
        N6_ind1=[-1;+1;-sizeIm1(1);+sizeIm1(1);-sizeIm1(1)*sizeIm1(2);+sizeIm1(1)*sizeIm1(2)];
        GDx1=GDx-min_x+pad1+1;
        GDy1=GDy-min_y+pad1+1;
        GDz1=GDz-min_z+pad1+1;
        GDind1=sub2ind(sizeIm1,GDx1,GDy1,GDz1);
        
        T1=inf(sizeIm1); T1(GDind1)=0;
        D1=inf(sizeIm1); D1(GDind1)=0;
        L1=zeros(sizeIm1); L1(GDind1)=1;
        KT1=zeros(sizeIm1); KT1(GDind1)=1;
                
        Xmin=max(min_x-1-pad1,0)+pad+1;
        Xmax=min(max_x-sizeIm(1)+pad1,0)+sizeIm(1)-pad;
        Ymin=max(min_y-1-pad1,0)+pad+1;
        Ymax=min(max_y-sizeIm(2)+pad1,0)+sizeIm(2)-pad;
        Zmin=max(min_z-1-pad1,0)+pad+1;
        Zmax=min(max_z-sizeIm(3)+pad1,0)+sizeIm(3)-pad;
        
        Xmin1=max(pad1-min_x+1,0)+pad+1;
        Xmax1=min(sizeIm(1)-max_x-pad1,0)+sizeIm1(1)-pad;
        Ymin1=max(pad1-min_y+1,0)+pad+1;
        Ymax1=min(sizeIm(2)-max_y-pad1,0)+sizeIm1(2)-pad;
        Zmin1=max(pad1-min_z+1,0)+pad+1;
        Zmax1=min(sizeIm(3)-max_z-pad1,0)+sizeIm1(3)-pad;
        
        Im1=zeros(sizeIm1);
        Im1(Xmin1:Xmax1,Ymin1:Ymax1,Zmin1:Zmax1)=Im(Xmin:Xmax,Ymin:Ymax,Zmin:Zmax);
        %Im1([1,2,end-1,end],:,:)=0; Im1(:,[1,2,end-1,end],:)=0; Im1(:,:,[1,2,end-1,end])=0;

        [T1,D1,~,~,~,~]=Fast_Marching(Im1,T1,D1,L1,KT1,N6_ind1,inf,Max_T);
        
        % merge the maps
        temp_ind=(T(Xmin:Xmax,Ymin:Ymax,Zmin:Zmax)<inf & L(Xmin:Xmax,Ymin:Ymax,Zmin:Zmax)==L(NewKnownPoint));
        T1=T1(Xmin1:Xmax1,Ymin1:Ymax1,Zmin1:Zmax1);
        D1=D1(Xmin1:Xmax1,Ymin1:Ymax1,Zmin1:Zmax1);
        T1(~temp_ind)=inf;
        D1(~temp_ind)=inf;
        T(Xmin:Xmax,Ymin:Ymax,Zmin:Zmax)=min(T(Xmin:Xmax,Ymin:Ymax,Zmin:Zmax),T1);
        D(Xmin:Xmax,Ymin:Ymax,Zmin:Zmax)=min(D(Xmin:Xmax,Ymin:Ymax,Zmin:Zmax),D1);
        
    elseif strcmp(exit_flag,'end')
        stop_condition=true;
    end
    
    figure(100)
    [NVx,NVy,NVz]=ind2sub(sizeIm,GDpath);
    plot3(NVy,NVx,NVz,'r.')
    drawnow
end

% Labels=unique(L(L>0));
% for i=1:length(Labels)
%     temp_ind1=(L==Labels(i));
%     temp_ind2=(T>0);
%     
%     L(temp_ind1 & temp_ind2)=0;
%     D(temp_ind1 & temp_ind2)=inf;
%     T(temp_ind1 & temp_ind2)=inf;
% end
L(T>0)=0;
%D(T>0)=inf;
%T(T>0)=inf;

disp('Fast Marching is complete.')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[EVx,EVy,EVz]=ind2sub(sizeIm,NewKnownPoint);

figure
imshow(max(L,[],3),[0 max(L(:))]), hold on
plot3(SVy,SVx,SVz,'r.')
plot3(EVy,EVx,EVz,'g.')

% figure
% Ltemp=L;
% Ltemp(KT==2)=0;
% imshow(max(Ltemp,[],3),[0 max(Ltemp(:))]), hold on
% plot3(SVy,SVx,SVz,'r.')
% plot3(EVy,EVx,EVz,'g.')
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
% this function performs Fast Marching untill stop condition is met
function [T,D,L,KT,exit_ind,exit_flag]=Fast_Marching(Im,T,D,L,KT,N6_ind,Max_dist,Max_time)
%KT=sparse(KT(:));
KnownPoints=find(KT==1);
% TrialPoints=find(KT==2);
% T_TrialPoints=T(TrialPoints);

NHood=[KnownPoints+N6_ind(1);KnownPoints+N6_ind(2);KnownPoints+N6_ind(3); ...
    KnownPoints+N6_ind(4);KnownPoints+N6_ind(5);KnownPoints+N6_ind(6)];
NHood=unique(NHood);
NHood(KT(NHood)==1)=[];

TrialPoints=NHood;
T_TrialPoints=T(TrialPoints);

stop_cond=true;
while stop_cond
    ind=[NHood+N6_ind(1),NHood+N6_ind(2),NHood+N6_ind(3),NHood+N6_ind(4),NHood+N6_ind(5),NHood+N6_ind(6)];
    
    temp_T=T(ind);
    temp_T(KT(ind)~=1)=inf;
    [T(NHood),D(NHood)]=BackSolver(Im(NHood),temp_T,D(ind));
    
    [~,temp_ind]=min(temp_T,[],2);
    
    temp=(1:length(NHood))'+length(NHood).*(temp_ind-1);
    L(NHood)=L(ind(temp));
    
    keep=NHood(KT(NHood)==0);
    KT(keep)=2;
    TrialPoints=[TrialPoints;keep];
    T_TrialPoints=[T_TrialPoints;T(keep)];
    
    [~,min_ind]=min(T_TrialPoints);
    NewKnownPoint=TrialPoints(min_ind);
    if ~isinf(NewKnownPoint)
        KT(NewKnownPoint)=1;
        T_TrialPoints(min_ind)=inf;
        NHood=NewKnownPoint+N6_ind;
        
        tempL=L(NHood(KT(NHood)==1));
        tempL=tempL(tempL>0);
        collision=(sum(tempL-tempL(1))~=0);
        NHood(KT(NHood)==1)=[];
        
        if D(NewKnownPoint)>=Max_dist
            stop_cond=false;
            exit_flag='dist'; % maximum distance reached
        end
        if T(NewKnownPoint)>=Max_time
            stop_cond=false;
            exit_flag='time'; % maximum time reached
        end
        if collision
            stop_cond=false; 
            exit_flag='collision'; % collision
        end
    else
        stop_cond=false;
        exit_flag='end'; % No place to propagate.  
    end
end
exit_ind=NewKnownPoint;

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
indd=((cumsum(Dpm_xyz.^2,2)-cumsum(Dpm_xyz,2).^2./N)<Cd);
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
function [T,D,L,KT,GDpath]=GradientDescent(T,D,L,KT,ind,N6_ind)

pot_ind=ind+N6_ind;
pot_ind(KT(pot_ind)~=1)=[];
Labels=unique(L(pot_ind));

GDpath=zeros(1,sum(size(T)));
GDpath(1)=ind;
count=1;
for i=1:length(Labels)
    current_ind=ind;
    while T(current_ind)~=0
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
            GDpath(count)=current_ind;
        end
    end
end
GDpath=GDpath(1:count);

% relabel
if length(Labels)>1
    L_new=max(L(:))+1;
    for i=1:length(Labels)
        L((L==Labels(i)))=L_new;
    end
else
    L_new=L(ind);
end
T(GDpath)=0;
D(GDpath)=0;
L(GDpath)=L_new;
KT(GDpath)=1;




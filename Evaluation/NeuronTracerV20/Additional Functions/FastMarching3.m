% This function solves of the Eikonal equation by using the Fast Marching algorithm 
% of Sethian. T is the arival time map. The boundary condition is set as 
% T(Start_X,Start_Y,Start_Z) = 0 (if start point is provided). 

function [Im,AMlbl,r]=FastMarching3(Orig)

disp('Fast Marching started.')
cl=class(Orig);

N_steps=1000; % number of re-initialization steps
Max_dist=10; % re-initialization is done if a fornt reaches Max_dist
Min_dist=5;
Max_time=inf;
seed_R=50; 
seed_thr=100;

%Orig=Multi_Scale_LoG(Orig,3,1,3);

[SVx,SVy,SVz,~]=Find_Seeds(Orig,seed_R,seed_thr);
%SVx=194; SVy=292; SVz=19;
%SVx=340; SVy=626; SVz=7;

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

AMlbl=zeros(length(StartVoxel),length(StartVoxel));
r_ind=StartVoxel;

N6_ind=[-1;+1;-sizeIm(1);+sizeIm(1);-sizeIm(1)*sizeIm(2);+sizeIm(1)*sizeIm(2)];

pad1=fix(1*Max_dist)+pad;

T=inf(sizeIm); T(StartVoxel)=0;
D=inf(sizeIm); D(StartVoxel)=0;
L=zeros(sizeIm); L(StartVoxel)=(1:length(StartVoxel));
KT=zeros(sizeIm); KT(StartVoxel)=1; % Known=1, Trial=2
KnownPoints=StartVoxel;

figure(100)
imshow(max(Im,[],3),[0 max(Im(:))]), hold on
[NVx,NVy,NVz]=ind2sub(sizeIm,StartVoxel);
plot3(NVy,NVx,NVz,'g.')
drawnow

stop_condition=false;
count=0;
Lcount=length(StartVoxel);
while stop_condition==false && count<=N_steps
    display([count,Lcount]);

    % do fast marching untill 1 of 4 conditions is encountered
    [T,D,L,KT,KnownPoints,NewKnownPoint,exit_flag]=Fast_Marching(Im,T,D,L,KT,KnownPoints,N6_ind,Max_dist,Max_time);
    
    Max_T=T(NewKnownPoint);
    if strcmp(exit_flag,'dist') || strcmp(exit_flag,'collision')
        % find the gradient descent path and relabel the regions
        [T,D,L,KT,AMlbl,r_ind,GDpath,Nregions]=GradientDescent(T,D,L,KT,AMlbl,r_ind,NewKnownPoint,N6_ind);
        count=count+1;
        Lcount=Lcount+1-Nregions; 
        
        [GDx,GDy,GDz]=ind2sub(sizeIm,GDpath);
        min_x=min(GDx); max_x=max(GDx);
        min_y=min(GDy); max_y=max(GDy);
        min_z=min(GDz); max_z=max(GDz);
                                       
        Xmin=max(min_x-1-pad1,0)+pad+1;
        Xmax=min(max_x-sizeIm(1)+pad1,0)+sizeIm(1)-pad;
        Ymin=max(min_y-1-pad1,0)+pad+1;
        Ymax=min(max_y-sizeIm(2)+pad1,0)+sizeIm(2)-pad;
        Zmin=max(min_z-1-pad1,0)+pad+1;
        Zmax=min(max_z-sizeIm(3)+pad1,0)+sizeIm(3)-pad;
        
        Im1=Im(Xmin:Xmax,Ymin:Ymax,Zmin:Zmax);
        Im1=padarray(Im1,[pad,pad,pad]);
        sizeIm1=size(Im1);
        N6_ind1=[-1;+1;-sizeIm1(1);+sizeIm1(1);-sizeIm1(1)*sizeIm1(2);+sizeIm1(1)*sizeIm1(2)];
        GDx1=GDx-Xmin+1+pad;
        GDy1=GDy-Ymin+1+pad;
        GDz1=GDz-Zmin+1+pad;
        GDind1=sub2ind(sizeIm1,GDx1,GDy1,GDz1);
        
        T1=inf(sizeIm1); T1(GDind1)=0;
        D1=inf(sizeIm1); D1(GDind1)=0;
        L1=zeros(sizeIm1); L1(GDind1)=1;
        KT1=zeros(sizeIm1); KT1(GDind1)=1;

        [T1,D1,~,~,~,~,~]=Fast_Marching(Im1,T1,D1,L1,KT1,GDind1,N6_ind1,inf,Max_T);
        
        % merge the maps
        temp_ind=(L(Xmin:Xmax,Ymin:Ymax,Zmin:Zmax)==L(NewKnownPoint));
        T1=T1(pad+1:end-pad,pad+1:end-pad,pad+1:end-pad);
        D1=D1(pad+1:end-pad,pad+1:end-pad,pad+1:end-pad);
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
    
%     rem_ind=(D>Min_dist & D<inf);
%     T(rem_ind)=inf;
%     D(rem_ind)=inf;
%     L(rem_ind)=0;
%     KT(rem_ind)=0;
%     KnownPoints=find(KT==1);
end

Im=zeros(size(Im),cl);
Im(T==0)=intmax(cl);
Im=Im(1+pad:end-pad,1+pad:end-pad,1+pad:end-pad);
[rx,ry,rz]=ind2sub(sizeIm,r_ind);
r=[rx,ry,rz]-pad;

rem=(sum(AMlbl,1)==0);
AMlbl(rem,:)=[];
AMlbl(:,rem)=[];
r(rem,:)=[];
AMlbl = LabelTreesAM(AMlbl);

disp('Fast Marching is complete.')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure
imshow(max(Orig,[],3),[0 max(Orig(:))])
hold on
PlotAM(AMlbl, r)
plot3(SVy-pad,SVx-pad,SVz-pad,'g.')

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
function [T,D,L,KT,KnownPoints,exit_ind,exit_flag]=Fast_Marching(Im,T,D,L,KT,KnownPoints,N6_ind,Max_dist,Max_time)

NHood=[KnownPoints+N6_ind(1);KnownPoints+N6_ind(2);KnownPoints+N6_ind(3); ...
    KnownPoints+N6_ind(4);KnownPoints+N6_ind(5);KnownPoints+N6_ind(6)];
NHood(KT(NHood)==1)=[];
NHood=unique(NHood);

TrialPoints=NHood;
T_TrialPoints=T(TrialPoints);

stop_cond=true;
while stop_cond
    ind=[NHood+N6_ind(1),NHood+N6_ind(2),NHood+N6_ind(3),NHood+N6_ind(4),NHood+N6_ind(5),NHood+N6_ind(6)];
    
    KT_1=(KT(ind)~=1);
    temp_T=T(ind);
    temp_T(KT_1)=inf;
    temp_D=D(ind);
    temp_D(KT_1)=inf;
    [T(NHood),D(NHood)]=BackSolver(Im(NHood),temp_T,temp_D);
    
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
        KnownPoints=[KnownPoints;NewKnownPoint];
        
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
function [T,D,L,KT,AMlbl,r_ind,GDpath,Nregions]=GradientDescent(T,D,L,KT,AMlbl,r_ind,ind,N6_ind)

lengthAM=size(AMlbl,1);
pot_ind=ind+N6_ind;
pot_ind(KT(pot_ind)~=1)=[];
Labels=unique(L(pot_ind));
Nregions=length(Labels);

GDpath=zeros(sum(size(T)),1);
GDpath(1)=ind;
AMlbl(lengthAM+1,lengthAM+1)=0;
count=1;
for i=1:Nregions
    current_ind=ind;
    while T(current_ind)~=0
        pot_ind=current_ind+N6_ind;
        pot_ind(L(pot_ind)~=Labels(i))=[];
        pot_ind(KT(pot_ind)~=1)=[];
        pot_ind(pot_ind==ind)=[];
        
        [~,temp]=min(T(pot_ind));
        current_ind=pot_ind(temp);
        if ~isempty(current_ind)
            if count==100
                error(num2str(count))
            end
            count=count+1;
            GDpath(count)=current_ind;
        end
    end
    count=count-1;
    Start=size(AMlbl,1)+1;
    End=lengthAM+count;
    if End>=Start
        AMlbl(End,End)=0;
        AMlbl(Start:End,Start:End)=diag(nan(1,End-Start),1)+diag(nan(1,End-Start),-1);
        AMlbl(lengthAM+1,Start)=NaN;
        AMlbl(Start,lengthAM+1)=NaN;
        temp=find(r_ind==current_ind);
        AMlbl(temp,End)=NaN;
        AMlbl(End,temp)=NaN;
    elseif End==Start-1
        temp=find(r_ind==current_ind);
        AMlbl(temp,lengthAM+1)=NaN;
        AMlbl(lengthAM+1,temp)=NaN;    
    end
end

GDpath=GDpath(1:count);
r_ind=[r_ind;GDpath];

% relabel
if Nregions>1
    L_new=max(L(:))+1;
    for i=1:Nregions
        L((L==Labels(i)))=L_new;
        AMlbl((AMlbl==Labels(i)))=L_new;
    end
else
    L_new=L(ind);
end

T(GDpath)=0;
D(GDpath)=0;
L(GDpath)=L_new;
KT(GDpath)=1;
AMlbl(isnan(AMlbl))=L_new;


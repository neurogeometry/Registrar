% This function solves of the Eikonal equation by using the Fast Marching algorithm 
% of Sethian. T is the arival time map. The boundary condition is set as 
% T(Start_X,Start_Y,Start_Z) = 0 (if start point is provided). 

function [Im,AMlbl,r]=FastMarching4(Orig)

disp('Fast Marching started.')
cl=class(Orig);

N_steps=1000; % number of re-initialization steps
Max_dist=15; % re-initialization is done if a fornt reaches Max_dist
Min_dist=inf;
Max_time=inf;
seed_R=30; 
seed_thr=50;
pad=2;
pad1=ceil(Max_dist)+pad;
background=0.0;
create_loops=false;

Orig1=Multi_Scale_LoG(Orig,2,1,3);
[SVx,SVy,SVz,~]=Find_Seeds(Orig1,seed_R,seed_thr);

sizeOrig=size(Orig);
if length(sizeOrig)==2
    sizeOrig=[sizeOrig,1];
end

Orig=double(Orig);
Orig=Orig-min(Orig(:));
Orig=(Orig./max(Orig(:)));
Orig(Orig<background)=background;

Im=zeros(sizeOrig+2*pad);
Im(1+pad:end-pad,1+pad:end-pad,1+pad:end-pad)=Orig;
sizeIm=size(Im);
SVx=SVx+pad;
SVy=SVy+pad;
SVz=SVz+pad;
StartVoxel=sub2ind_AS(sizeIm,SVx,SVy,SVz);

AMlbl=sparse(length(StartVoxel),length(StartVoxel));
r_ind=StartVoxel;

N6_ind=[-1;+1;-sizeIm(1);+sizeIm(1);-sizeIm(1)*sizeIm(2);+sizeIm(1)*sizeIm(2)];

T=inf(sizeIm); T(StartVoxel)=0;
D=inf(sizeIm); D(StartVoxel)=0;
%I=inf(sizeIm); I(StartVoxel)=Im(StartVoxel);
L=zeros(sizeIm); L(StartVoxel)=(1:length(StartVoxel));
AML=sparse(diag(ones(1,1000)));
KT=zeros(sizeIm); KT(StartVoxel)=1; % Known=1, Trial=2
KnownPoints=StartVoxel;

figure(100)
imshow(max(Im,[],3),[0 max(Im(:))]), hold on
[NVx,NVy,NVz]=ind2sub(sizeIm,StartVoxel);
plot3(NVy,NVx,NVz,'g.')
drawnow

stop_condition=false;
Q=zeros(1,N_steps);
Ntrace=zeros(1,N_steps);
Nknown=zeros(1,N_steps);
Itrace=zeros(1,N_steps);
Iknown=zeros(1,N_steps);
% Q1=zeros(1,N_steps);
% Ntrace1=zeros(1,N_steps);
% Nknown1=zeros(1,N_steps);
% Itrace1=zeros(1,N_steps);
% Iknown1=zeros(1,N_steps);
av_window=20;
count=0;
Lcount=length(StartVoxel);
while stop_condition==false && count<=N_steps
    display([count,Lcount]);

    % do fast marching untill 1 of 4 conditions is encountered
    [T,D,L,AML,KT,KnownPoints,NewKnownPoint,exit_flag]=Fast_Marching(Im,T,D,L,AML,KT,KnownPoints,N6_ind,Max_dist,Max_time);
    
    Max_T=T(NewKnownPoint);
    if strcmp(exit_flag,'dist') || strcmp(exit_flag,'collision')
        % find the gradient descent path and relabel the regions
        [T,D,L,AML,KT,AMlbl,r_ind,GDpath,Nregions]=GradientDescent(T,D,L,AML,KT,AMlbl,r_ind,NewKnownPoint,N6_ind,Max_dist,create_loops);
        if length(AMlbl)~=length(r_ind)
            1
        end
        count=count+1;
        Lcount=Lcount+1-Nregions;
        
        if ~isempty(GDpath)
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
            %I1=inf(sizeIm1); I1(GDind1)=Im1(GDind1);
            L1=zeros(sizeIm1); L1(GDind1)=L(GDpath(1));
            AML1=sparse(diag(ones(1,1000)));
            KT1=zeros(sizeIm1); KT1(GDind1)=1;
            
            [T1,D1,L1,~,~,~,~,~]=Fast_Marching(Im1,T1,D1,L1,AML1,KT1,GDind1,N6_ind1,inf,Max_T);
            
            %         Ntrace1(count)=sum(T1(:)==0);
            %         Nknown1(count)=sum(KT1(:)==1);
            %         Itrace1(count)=sum(Im1(T1(:)==0));
            %         Iknown1(count)=sum(Im1(KT1(:)==1));
            % Nknown1./Itrace1, Nknown1./Ntrace1
            
            % merge the maps
            T1=T1(pad+1:end-pad,pad+1:end-pad,pad+1:end-pad);
            D1=D1(pad+1:end-pad,pad+1:end-pad,pad+1:end-pad);
            L1=L1(pad+1:end-pad,pad+1:end-pad,pad+1:end-pad);
            %I1=I1(pad+1:end-pad,pad+1:end-pad,pad+1:end-pad);
            
            %         pot_ind=NewKnownPoint+N6_ind;
            %         pot_ind(KT(pot_ind)~=1)=[];
            %         pot_Labels=unique(L(pot_ind));
            pot_Labels=find(AML(L(NewKnownPoint),:));
            
            % overwrite existing maps if:   L==pot_Labels & T1<=T
            % or:                           L~=pot_Labels & (KT~=1 & T1<=T)
            temp_ind=zeros(size(T1));
            for j=1:length(pot_Labels)
                temp_ind=(temp_ind | (L(Xmin:Xmax,Ymin:Ymax,Zmin:Zmax)==pot_Labels(j)));
            end
            
            T1(~temp_ind & KT(Xmin:Xmax,Ymin:Ymax,Zmin:Zmax)==1)=inf;
            D1(~temp_ind & KT(Xmin:Xmax,Ymin:Ymax,Zmin:Zmax)==1)=inf;
            L1(~temp_ind & KT(Xmin:Xmax,Ymin:Ymax,Zmin:Zmax)==1)=inf;
            %I1(~temp_ind)=inf;
            
            %         temp_ind=(T1>T(Xmin:Xmax,Ymin:Ymax,Zmin:Zmax));
            %         T1(temp_ind)=inf;
            %         D1(temp_ind)=inf;
            %         L1(temp_ind)=inf;
            %         T(Xmin:Xmax,Ymin:Ymax,Zmin:Zmax)=min(T(Xmin:Xmax,Ymin:Ymax,Zmin:Zmax),T1);
            %         D(Xmin:Xmax,Ymin:Ymax,Zmin:Zmax)=min(D(Xmin:Xmax,Ymin:Ymax,Zmin:Zmax),D1);
            %         L(Xmin:Xmax,Ymin:Ymax,Zmin:Zmax)=min(L(Xmin:Xmax,Ymin:Ymax,Zmin:Zmax),L1);
            %I(Xmin:Xmax,Ymin:Ymax,Zmin:Zmax)=min(I(Xmin:Xmax,Ymin:Ymax,Zmin:Zmax),I1);
            
            %         Ntrace1(count)=sum(T1(:)==0);
            %         Nknown1(count)=sum(KT1(:)==1);
            %         Itrace1(count)=sum(Im1(T1(:)==0));
            %         Iknown1(count)=sum(Im1(KT1(:)==1));
            % Nknown1./Itrace1, Nknown1./Ntrace1
            
            temp_ind=(T1<=T(Xmin:Xmax,Ymin:Ymax,Zmin:Zmax));
            tempT=T(Xmin:Xmax,Ymin:Ymax,Zmin:Zmax); tempT(temp_ind)=T1(temp_ind); T(Xmin:Xmax,Ymin:Ymax,Zmin:Zmax)=tempT;
            tempD=D(Xmin:Xmax,Ymin:Ymax,Zmin:Zmax); tempD(temp_ind)=D1(temp_ind); D(Xmin:Xmax,Ymin:Ymax,Zmin:Zmax)=tempD;
            tempL=L(Xmin:Xmax,Ymin:Ymax,Zmin:Zmax); tempL(temp_ind)=L1(temp_ind); L(Xmin:Xmax,Ymin:Ymax,Zmin:Zmax)=tempL;
        end
    elseif strcmp(exit_flag,'end')
        stop_condition=true;
    end
    
    figure(100)
    [NVx,NVy,NVz]=ind2sub(sizeIm,GDpath);
    plot3(NVy,NVx,NVz,'r.')
    drawnow
    
    Ntrace(count)=sum(T(:)==0);
    Nknown(count)=sum(KT(:)==1);
    Itrace(count)=sum(Im(T(:)==0));
    Iknown(count)=sum(Im(KT(:)==1));
    % Nknown./Itrace, Nknown./Ntrace 
    
    %Q1(count)=Nknown1(count)/Itrace1(count);
    if count>=2
        %Q(count)=(Nknown(count)/Ntrace(count)/pi)^0.5-(Nknown(count-1)/Ntrace(count-1)/pi)^0.5;
        Q(count)=Nknown(count)/Itrace(count)-Nknown(count-1)/Itrace(count-1);
        if count>av_window
            Qmean=nanmean(Q(count-av_window:count-1));
            Qstd=nanstd(Q(count-av_window:count-1));
            %Q1mean=nanmean(Q1(count-av_window:count-1));
            %Q1std=nanstd(Q1(count-av_window:count-1));
        else
            Qmean=nanmean(Q(1:count-1));
            Qstd=nanstd(Q(1:count-1));
            %Q1mean=nanmean(Q1(1:count-1));
            %Q1std=nanstd(Q1(1:count-1));
        end
        %[(Q(count)-Qmean)/Qstd,(Q1(count)-Q1mean)/Q1std]
        if count>av_window && Q(count)>Qmean+5*Qstd %|| Q1(count)>Q1mean+5*Q1std)
            reply = input('Do you want to continue? ', 's');
            if strcmp(reply,'y')
                Q(count)=NaN;
                %Q1(count)=NaN;
            else
                stop_condition=true;
            end
        end
    end
    
    if Min_dist<=Max_dist
        rem_ind=(D>Min_dist & D<inf);
        T(rem_ind)=inf;
        D(rem_ind)=inf;
        L(rem_ind)=0;
        KT(rem_ind)=0;
        KnownPoints=find(KT==1);
    end
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
AMlbl = LabelBranchesAM(AMlbl);
[AMlbl, r, ~, ~]=Reduce_Small_Loops(AMlbl,r,Max_dist);
AMlbl = LabelTreesAM(AMlbl);

disp('Fast Marching is complete.')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure(300)
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
function [T,D,L,AML,KT,KnownPoints,exit_ind,exit_flag]=Fast_Marching(Im,T,D,L,AML,KT,KnownPoints,N6_ind,Max_dist,Max_time)

NHood=[KnownPoints+N6_ind(1);KnownPoints+N6_ind(2);KnownPoints+N6_ind(3); ...
    KnownPoints+N6_ind(4);KnownPoints+N6_ind(5);KnownPoints+N6_ind(6)];
NHood(KT(NHood)==1)=[];
NHood=unique(NHood);

TrialPoints=NHood;
T_TrialPoints=T(TrialPoints);

stop_cond=true;
if isempty(NHood)
    stop_cond=false;
    exit_flag='end';
end

while stop_cond
    ind=[NHood+N6_ind(1),NHood+N6_ind(2),NHood+N6_ind(3),NHood+N6_ind(4),NHood+N6_ind(5),NHood+N6_ind(6)];

    KT_1=(KT(ind)~=1);
    temp_T=T(ind);
    temp_T(KT_1)=inf;
    temp_D=D(ind);
    temp_D(KT_1)=inf;
    %temp_I=I(ind);
    %temp_I(KT_1)=inf;
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
    if ~isinf(T(NewKnownPoint)) && ~isempty(NewKnownPoint)
        KT(NewKnownPoint)=1;
        T_TrialPoints(min_ind)=inf;
        NHood=NewKnownPoint+N6_ind;
        KnownPoints=[KnownPoints;NewKnownPoint];
        
        tempL=L(NHood(KT(NHood)==1));
        tempL=unique(tempL(tempL>0));
        %collision=(sum(tempL-L(NewKnownPoint))~=0);
        collision=(sum(AML(L(NewKnownPoint),tempL))~=length(tempL));

        if collision==1
            if sum(AML(L(NewKnownPoint),tempL))==length(tempL)
                collision=0;
            else
                AML(L(NewKnownPoint),tempL)=1;
                AML(tempL,L(NewKnownPoint))=1;
            end
        end
        NHood(KT(NHood)==1)=[];
        
        if T(NewKnownPoint)>=Max_time
            stop_cond=false;
            exit_flag='time'; % maximum time reached
        end
        if D(NewKnownPoint)>=Max_dist
            stop_cond=false;
            exit_flag='dist'; % maximum distance reached
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
% this function finds arrival time, T, distance, D, and integrated intensity, I, for the trial points
function [T,D]=BackSolver(F,Tpm_xyz,Dpm_xyz)

delta=1;
N=ones(size(Tpm_xyz,1),1)*(1:size(Tpm_xyz,2));
ct=(delta./F).^2;
cd=(delta).^2.*ones(size(F));
%ci=(delta.*F).^2;
Ct=ct*ones(1,size(Tpm_xyz,2));
Cd=cd*ones(1,size(Dpm_xyz,2));
%Ci=ci*ones(1,size(Ipm_xyz,2));
Tpm_xyz=sort(Tpm_xyz,2);
Dpm_xyz=sort(Dpm_xyz,2);
%Ipm_xyz=sort(Ipm_xyz,2);
indt=((cumsum(Tpm_xyz.^2,2)-cumsum(Tpm_xyz,2).^2./N)<=Ct);
indd=((cumsum(Dpm_xyz.^2,2)-cumsum(Dpm_xyz,2).^2./N)<Cd);
%indi=((cumsum(Ipm_xyz.^2,2)-cumsum(Ipm_xyz,2).^2./N)<=Ci);
Tpm_xyz(~indt)=0;
Dpm_xyz(~indd)=0;
%Ipm_xyz(~indi)=0;
nt=sum(indt,2);
nd=sum(indd,2);
%ni=sum(indi,2);
Tmean=sum(Tpm_xyz,2)./nt;
Dmean=sum(Dpm_xyz,2)./nd;
%Imean=sum(Ipm_xyz,2)./ni;
T2mean=sum(Tpm_xyz.^2,2)./nt;
D2mean=sum(Dpm_xyz.^2,2)./nd;
%I2mean=sum(Ipm_xyz.^2,2)./ni;
T=Tmean+(ct./nt-(T2mean-Tmean.^2)).^0.5;
D=Dmean+(cd./nd-(D2mean-Dmean.^2)).^0.5;
%I=Imean+(ci./ni-(I2mean-Imean.^2)).^0.5;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function performs gradient descent from ind (a single point) to
% T=0. T map of the region is reset to the gradient descent path, Tnew.
% If there are multiple labels neighboring ind, the gradient descent will
% go into all these labeld regions. The result will be merged and
% relabeled.
function [T,D,L,AML,KT,AMlbl,r_ind,GDpath,Nregions]=GradientDescent(T,D,L,AML,KT,AMlbl,r_ind,ind,N6_ind,Max_dist,create_loops)

lengthAM=size(AMlbl,1);
pot_ind=ind+N6_ind;
pot_ind(KT(pot_ind)~=1)=[];
Labels=unique(L(pot_ind));
% Ttemp=T(3:end-2,3:end-2,3:end-2);
% Ltemp=L(3:end-2,3:end-2,3:end-2);
% for i=1:length(Labels)
%     if isempty(find(Ttemp(Ltemp==Labels(i))==0,1,'first'))
%         Labels(i)=NaN;
%     end
% end
% Labels(isnan(Labels))=[];
Nregions=length(Labels);

GDpath=zeros(sum(size(T)),1);
GDpath(1)=ind;
AMlbl(lengthAM+1,lengthAM+1)=0;
count=1;
for i=1:Nregions
    current_ind=ind;
    count1=1;
    while ~isempty(current_ind)&& T(current_ind)~=0 && count1<=2*Max_dist
        pot_ind=current_ind+N6_ind;
        pot_ind(L(pot_ind)~=Labels(i))=[];
        pot_ind(KT(pot_ind)~=1)=[];
        pot_ind(pot_ind==ind)=[];
        
        [~,temp]=min(T(pot_ind));
        current_ind=pot_ind(temp);
        count1=count1+1;
        if ~isempty(current_ind)
            count=count+1;
            GDpath(count)=current_ind;
        end
    end
    
    current_ind=GDpath(count);
    count=count-1;
    Start=size(AMlbl,1)+1;
    End=lengthAM+count;
    if count1<=2*Max_dist   
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
    else
        disp('Gradient descent failed to find a T=0 voxel.')
        if End>=Start
            AMlbl(End,End)=0;
        elseif End==Start-1
            temp=find(r_ind==current_ind);
            AMlbl(temp,lengthAM+1)=NaN;
            AMlbl(lengthAM+1,temp)=NaN;
        end
    end
end

GDpath=GDpath(1:count);
GDpath(isnan(GDpath))=[];
r_ind=[r_ind;GDpath];

% relabel
if create_loops
    L_new=max(L(:))+1;
    AML(L_new,Labels)=1;
    AML(Labels,L_new)=1;
    AML(L_new,L_new)=1;
else
    if Nregions>1
        L_new=max(L(:))+1;
        for i=1:Nregions
            L((L==Labels(i)))=L_new;
            AMlbl((AMlbl==Labels(i)))=L_new;
        end
    else
        L_new=L(ind);
    end
end

T(GDpath)=0;
D(GDpath)=0;
%L(GDpath)=L_new;
KT(GDpath)=1;
AMlbl(isnan(AMlbl))=L_new;


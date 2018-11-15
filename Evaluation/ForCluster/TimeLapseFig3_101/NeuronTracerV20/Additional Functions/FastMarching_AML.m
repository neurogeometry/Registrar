% This function solves of the Eikonal equation by using the Fast Marching algorithm
% of Sethian. T and D are the arival time and distance maps.
% N_steps is the maximum number of re-initialization steps
% Max_Known_Dist is the distance at which re-initialization is performed
% Min_Known_Dist is the distance to which map is truncated after re-initialization
% Contrast is the minimum ratio of added branch and trace intensities  
% SVr contains positions of the seeds
% unisotropy is the wave speed unisotropy in a uniform intensity image  
% This version uses AM for labels to determine collisions.
% GD problem

function [AMlbl,r]=FastMarching_AML(Orig,SVr,N_steps,Max_Known_Dist,Contrast,R_min,R_step,R_max,unisotropy,varargin)

if length(varargin)==1
    myClass=varargin{1};
else
    myClass.getFlag=NaN;
end

create_loops=true;
pad=2;
Min_Known_Dist=inf;

unisotropy(unisotropy<0)=0;
if sum(unisotropy)==0
    unisotropy=[1,1,1];
end
unisotropy=unisotropy./sum(unisotropy.^2)^0.5;
h1=1./[unisotropy(1),unisotropy(1),unisotropy(2),unisotropy(2),unisotropy(3),unisotropy(3)]';
h2=[unisotropy(1),unisotropy(1),unisotropy(2),unisotropy(2),unisotropy(3),unisotropy(3)].^2;

disp(['   # steps  ', '# contrast'])

sizeOrig=size(Orig);
if length(sizeOrig)==2 
    sizeOrig=[sizeOrig,1];
end

Orig=double(Orig);
Orig=Orig-min(Orig(:));
Orig=(Orig./max(Orig(:)));
background=mean(Orig(:));
Orig(Orig<background)=background;
%Orig=Orig+rand(sizeOrig).*background;

Im=zeros(sizeOrig+2*pad);
Im(1+pad:end-pad,1+pad:end-pad,1+pad:end-pad)=Orig;
clear Orig
sizeIm=sizeOrig+2*pad;
SVr=SVr+pad;
StartVoxel=sub2ind_AS(sizeIm,SVr(:,1),SVr(:,2),SVr(:,3));

AMlbl=sparse(length(StartVoxel),length(StartVoxel));
r_ind=StartVoxel;
[R I] = Find_R_fast(Im,SVr,R_min,max(0.5,(R_max-R_min)/10),R_max);

N6_ind=[-1;+1;-sizeIm(1);+sizeIm(1);-sizeIm(1)*sizeIm(2);+sizeIm(1)*sizeIm(2)];

T=inf(sizeIm); T(StartVoxel)=0;
D=inf(sizeIm); D(StartVoxel)=0;

L=zeros(sizeIm,'uint16'); L(StartVoxel)=(1:length(StartVoxel));
Lmax=length(StartVoxel);
AML=sparse(diag(ones(1,2*length(StartVoxel))));
KT=zeros(sizeIm,'uint8'); KT(StartVoxel)=1; % Known=1, Trial=2
KnownPoints=StartVoxel;

if isnan(myClass.getFlag)
    figure(100)
    imshow(max(Im,[],3),[0 max(Im(:))]), hold on
    [NVx,NVy,NVz]=ind2sub(sizeIm,StartVoxel);
    plot3(NVy,NVx,NVz,'r.','MarkerSize',20)
    drawnow
end

stop_condition=false;
count=0;
CC=NaN;
pad1=ceil(Max_Known_Dist)+pad;
while stop_condition==false   
    % flag values:
    % 1=Cancel
    % 2=Stop
    % 3=Apply
    % 4=Run/Continue
    
    while myClass.getFlag()==2
        pause(1);
    end
        
    if myClass.getFlag==3
        [rx,ry,rz]=ind2sub(sizeIm,r_ind);
        r=[rx,ry,rz]-pad;
        
        rem=(sum(AMlbl,2)==0 | isnan(r_ind));
        AMlbl(rem,:)=[];
        AMlbl(:,rem)=[];
        r(rem,:)=[];
        return
    elseif myClass.getFlag==1
        AMlbl=[];
        r=[];
        return
    elseif myClass.getFlag==4
        N_steps=myClass.getN_steps();
        Max_Known_Dist=myClass.getMax_Known_Dist();
        pad1=ceil(Max_Known_Dist)+pad;
        Contrast=myClass.getContrast();
        
        SVr_new=myClass.getSVr()+pad;
        if ~isempty(SVr_new)
            StartVoxel_new=sub2ind_AS(sizeIm,SVr_new(:,1),SVr_new(:,2),SVr_new(:,3));
            AMlbl(size(AMlbl,1)+length(StartVoxel_new),size(AMlbl,2)+length(StartVoxel_new))=0;
            r_ind=[r_ind;StartVoxel_new];
            [R_new I_new] = Find_R_fast(Im,SVr_new,R_min,max(0.5,(R_max-R_min)/10),R_max);
            R=[R;R_new];
            I=[I,I_new];
            
            T(StartVoxel_new)=0;
            D(StartVoxel_new)=0;
            L(StartVoxel_new)=Lmax+(1:length(StartVoxel_new));
            Lmax=Lmax+length(StartVoxel_new);
            AML(end+1:end+2*length(StartVoxel_new),end+1:end+2*length(StartVoxel_new))=sparse(diag(ones(1,2*length(StartVoxel_new))));
            KT(StartVoxel_new)=1;
            KnownPoints=[KnownPoints;StartVoxel_new];
        end
    end
       
    % do fast marching untill 1 of 4 conditions is encountered
    [T,D,L,AML,KT,KnownPoints,NewKnownPoint,exit_flag]=Fast_Marching_Step(Im,T,D,L,AML,KT,KnownPoints,N6_ind,Max_Known_Dist,inf,h2);
    
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     for ii=1:Lmax
%         temp=T(L==ii);
%         if ~isempty(temp) && sum(temp==0)==0
%             1
%         end
%     end
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    if strcmp(exit_flag,'dist') || strcmp(exit_flag,'collision')    
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % find the gradient descent path
        [AMlbl,r_ind,R,I,GDpath,~,I_GDpath,Labels]=GradientDescent(Im,T,L,KT,AMlbl,r_ind,R,R_min,R_step,R_max,I,NewKnownPoint,exit_flag,N6_ind,Max_Known_Dist,h1);
        
        if ~isempty(GDpath)
            CC=NaN;
            if strcmp(exit_flag,'dist') % add 'collision'
                CC=mean(I_GDpath(1:ceil(length(GDpath)/2)))/mean(I);
%               start=1;
%               endd=max(1,length(I_GDpath)-ceil(R(r_ind==GDpath(end))));
%               CC=mean(I_GDpath(start:endd))/mean(I);
%           elseif strcmp(exit_flag,'collision') && ~isempty(GDpath)
%               start=min(length(I_GDpath),ceil(R(r_ind==GDpath(end))));
%               endd=max(1,length(I_GDpath)-ceil(R(r_ind==GDpath(end))));
%               if start<endd
%                   CC=mean(I_GDpath(start:endd))/mean(I);
%               end
            end
        
            count=count+1;
            disp([count,CC])
                
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % relabel the regions
            T(GDpath)=0;
            D(GDpath)=0;
            KT(GDpath)=1;
            
            if create_loops 
                if strcmp(exit_flag,'dist')
                    L_new=Lmax+1;
                    Lmax=L_new;
                    AML(L_new,Labels)=1;
                    AML(Labels,L_new)=1;
                elseif strcmp(exit_flag,'collision')
                    L_new=Lmax+1;
                end
            elseif ~create_loops
                if length(Labels)>1
                    L_new=Lmax+1;
                    Lmax=L_new;
                    for i=1:length(Labels)
                        L((L==Labels(i)))=L_new;
                        AMlbl((AMlbl==Labels(i)))=L_new;
                    end
                else
                    L_new=L(NewKnownPoint);
                end
                L(GDpath)=L_new;
            end
            AMlbl(isnan(AMlbl))=L_new;
            
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
            
            T1old=T(Xmin:Xmax,Ymin:Ymax,Zmin:Zmax);
            T1old=padarray(T1old,[pad,pad,pad],inf);
            GDind1old=(T1old==0);
            
            T1=inf(sizeIm1); T1(GDind1old)=0; T1(GDind1)=0;
            D1=inf(sizeIm1); D1(GDind1old)=0; D1(GDind1)=0;
            L1=zeros(sizeIm1,'uint16'); L1(GDind1old)=1; L1(GDind1)=1;
            AML1=1;%sparse(diag(ones(1,2*length(AML))));
            KT1=zeros(sizeIm1,'uint8'); KT1(GDind1old)=1; KT1(GDind1)=1;
            
            [T1,D1,~,~,~,~,~,~,]=Fast_Marching_Step(Im1,T1,D1,L1,AML1,KT1,GDind1,N6_ind1,inf,max(T(KnownPoints)),h2);
                        
            % merge the maps
            T1=T1(pad+1:end-pad,pad+1:end-pad,pad+1:end-pad);
            D1=D1(pad+1:end-pad,pad+1:end-pad,pad+1:end-pad);
            %L1=L1(pad+1:end-pad,pad+1:end-pad,pad+1:end-pad);
            
            % overwrite existing maps if:   L==pot_Labels & T1<=T
            % or:                           L~=pot_Labels & (KT~=1 & T1<=T)

            pot_Labels=find(AML(L(NewKnownPoint),:));            
            temp_ind=zeros(size(T1));
            for j=1:length(pot_Labels)
                temp_ind=(temp_ind | (L(Xmin:Xmax,Ymin:Ymax,Zmin:Zmax)==pot_Labels(j)));
            end          
            T1(~temp_ind & KT(Xmin:Xmax,Ymin:Ymax,Zmin:Zmax)>-1)=inf;
            D1(~temp_ind & KT(Xmin:Xmax,Ymin:Ymax,Zmin:Zmax)>-1)=inf;
            %L1(~temp_ind & KT(Xmin:Xmax,Ymin:Ymax,Zmin:Zmax)==1)=inf;

            temp_indt=(T1<T(Xmin:Xmax,Ymin:Ymax,Zmin:Zmax) & KT(Xmin:Xmax,Ymin:Ymax,Zmin:Zmax)>0);
            tempT=T(Xmin:Xmax,Ymin:Ymax,Zmin:Zmax); tempT(temp_indt)=T1(temp_indt); T(Xmin:Xmax,Ymin:Ymax,Zmin:Zmax)=tempT;
            temp_indd=(D1<D(Xmin:Xmax,Ymin:Ymax,Zmin:Zmax) & KT(Xmin:Xmax,Ymin:Ymax,Zmin:Zmax)>0);
            tempD=D(Xmin:Xmax,Ymin:Ymax,Zmin:Zmax); tempD(temp_indd)=D1(temp_indd); D(Xmin:Xmax,Ymin:Ymax,Zmin:Zmax)=tempD;
            if strcmp(exit_flag,'dist')
                tempL=L(Xmin:Xmax,Ymin:Ymax,Zmin:Zmax); tempL(temp_indt)=L(GDpath(1)); L(Xmin:Xmax,Ymin:Ymax,Zmin:Zmax)=tempL;
            end
        end
        
        if Min_Known_Dist<=Max_Known_Dist
            rem_ind=(D>Min_Known_Dist & D<inf);
            T(rem_ind)=inf;
            D(rem_ind)=inf;
            L(rem_ind)=0;
            KT(rem_ind)=0;
            KnownPoints=find(KT==1);
        end
    end
        
    [NVx,NVy,NVz]=ind2sub(sizeIm,GDpath(~isnan(GDpath)));
    if isnan(myClass.getFlag)
        figure(100)
        plot3(NVy,NVx,NVz,'g.','MarkerSize',6)
        drawnow
    else
        ry=NVy-pad;
        rx=NVx-pad;
        rz=NVz-pad;
        myClass.updatePlot(rx,ry,rz);
    end
    
    if count>=N_steps
        disp('Maximum number of iteration steps reached. You may increase the number of iteration steps and continue Fast Marching.')
        if isnan(myClass.getFlag)
            stop_condition=true;
        else
            myClass.setFlag(2)    
        end
    end
    
    if CC<Contrast
        if isnan(myClass.getFlag)
            reply = input('Contrast of the last added branch is low. Do you want to continue? ', 's');
            if strcmp(reply,'n')
                stop_condition=true;
            elseif strcmp(reply,'y')
                stop_condition=false;
            end
        else
            myClass.setFlag(2);
            disp('Contrast of the last added branch is low. Click on "Continue" if you would like to resume Fast Marching.')
        end
    end
    
    if strcmp(exit_flag,'end')
        disp('Fast Marching has no place to propagate. You may add new seed points and continue Fast Marching.')
        if isnan(myClass.getFlag)
            stop_condition=true;
        else
            myClass.setFlag(2)
        end
    end
end

[rx,ry,rz]=ind2sub(sizeIm,r_ind);
r=[rx,ry,rz]-pad;

rem=(sum(AMlbl,2)==0 | isnan(r_ind));
AMlbl(rem,:)=[];
AMlbl(:,rem)=[];
r(rem,:)=[];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% this function performs Fast Marching untill a stop condition is met
function [T,D,L,AML,KT,KnownPoints,exit_ind,exit_flag]=Fast_Marching_Step(Im,T,D,L,AML,KT,KnownPoints,N6_ind,Max_Known_Dist,Max_Known_Time,h2)

NHood=[KnownPoints+N6_ind(1);KnownPoints+N6_ind(2);KnownPoints+N6_ind(3); ...
    KnownPoints+N6_ind(4);KnownPoints+N6_ind(5);KnownPoints+N6_ind(6)];
NHood(KT(NHood)==1)=[];
NHood(Im(NHood)==0)=[]; %!!!!!!!!!!!
NHood=unique(NHood);

TrialPoints=NHood(~isinf(T(NHood)));
T_TrialPoints=T(TrialPoints);

stop_cond=true;
if isempty(NHood)
    stop_cond=false;
    exit_flag='end';
end

while stop_cond
    if ~isempty(NHood)
        ind=[NHood+N6_ind(1),NHood+N6_ind(2),NHood+N6_ind(3),NHood+N6_ind(4),NHood+N6_ind(5),NHood+N6_ind(6)];
        KT_1=(KT(ind)~=1);
        
        Tpm_xyz=T(ind);
        Tpm_xyz(KT_1)=inf;
        [Tpm_xyz,IX]=sort(Tpm_xyz,2);
        
        Dpm_xyz=D(ind);
        Dpm_xyz(KT_1)=inf;
        Dpm_xyz=sort(Dpm_xyz,2);
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % this function finds arrival time, T, and distance, D, for the trial points
        
        H=h2(IX);
        H_cum=cumsum(H,2);
        ct=1./Im(NHood).^2;
        Tpm_xyz_cum=cumsum(Tpm_xyz.*H,2);
        Tpm_xyz2_cum=cumsum(Tpm_xyz.^2.*H,2);
        nt=sum(((Tpm_xyz2_cum-Tpm_xyz_cum.^2./H_cum)<=ct*ones(1,size(Tpm_xyz,2))),2);
        if sum(h2==0)>0
            nt(nt==0)=1;
        end
        ind_nt=(1:size(Tpm_xyz,1))'+(nt-1).*size(Tpm_xyz,1);
        ntH=H_cum(ind_nt);
        temp=Tpm_xyz_cum(ind_nt);
        T(NHood)=(temp+(ct.*ntH-(Tpm_xyz2_cum(ind_nt).*ntH-temp.^2)).^0.5)./ntH;
        
        N=ones(size(Tpm_xyz,1),1)*(1:size(Tpm_xyz,2));
        Dpm_xyz_cum=cumsum(Dpm_xyz,2);
        Dpm_xyz2_cum=cumsum(Dpm_xyz.^2,2);
        nd=sum(((Dpm_xyz2_cum-Dpm_xyz_cum.^2./N)<=ones(size(Dpm_xyz))),2);
        ind_nd=(1:size(Dpm_xyz,1))'+(nd-1).*size(Dpm_xyz,1);
        temp=Dpm_xyz_cum(ind_nd);
        D(NHood)=(temp+(nd-(Dpm_xyz2_cum(ind_nd).*nd-temp.^2)).^0.5)./nd;
        
        temp=(1:length(NHood))'+length(NHood).*(IX(:,1)-1);
        L(NHood)=L(ind(temp));
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        keep=NHood(KT(NHood)==0 & ~isinf(T(NHood)));
        KT(keep)=2;
        TrialPoints=[TrialPoints;keep];
        T_TrialPoints=[T_TrialPoints;T(keep)];
    end
    
    [~,min_ind]=min(T_TrialPoints);
    NewKnownPoint=TrialPoints(min_ind);
    
    if ~isempty(NewKnownPoint)
        KT(NewKnownPoint)=1;
        T_TrialPoints(min_ind)=inf;
        NHood=NewKnownPoint+N6_ind;
        KnownPoints=[KnownPoints;NewKnownPoint];
        
        tempL=L(NHood(KT(NHood)==1));
        collision=(sum(AML(L(NewKnownPoint),tempL))~=length(tempL));
        
        NHood(KT(NHood)==1)=[];
        NHood(Im(NHood)==0)=[]; %!!!!!!!!!!!
        
        if T(NewKnownPoint)>=Max_Known_Time
            stop_cond=false;
            exit_flag='time'; % maximum time reached
        end
        if D(NewKnownPoint)>=Max_Known_Dist
            stop_cond=false;
            exit_flag='dist'; % maximum distance reached
        end
        if collision
            stop_cond=false;
            exit_flag='collision';
            AML(L(NewKnownPoint),tempL)=1;
            AML(tempL,L(NewKnownPoint))=1;
            %AML(tempL,tempL)=1;
        end
    else
        stop_cond=false;
        exit_flag='end'; % No place to propagate.
    end
end
exit_ind=NewKnownPoint;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function performs gradient descent from ind (a single point) to
% T=0. T map of the region is reset to the gradient descent path, Tnew.
% If there are multiple labels neighboring ind, the gradient descent will
% go into all these labeld regions. The result will be merged and
% relabeled.
function [AMlbl,r_ind,R,I,GDpath,R_GDpath,I_GDpath,Labels]=GradientDescent(Im,T,L,KT,AMlbl,r_ind,R,R_min,R_step,R_max,I,ind,exit_flag,N6_ind,Max_Known_Dist,h1)

lengthAM=size(AMlbl,1);
pot_ind=ind+N6_ind;
pot_ind(KT(pot_ind)~=1)=[];
if strcmp(exit_flag,'collision')
    Labels=double(unique(L(pot_ind)));
else
    Labels=L(ind);
end

GDpath=zeros(sum(size(T)),1);
GDpath(1)=ind;
AMlbl(lengthAM+1,lengthAM+1)=0;
count=1;
for i=1:length(Labels)
    current_ind=ind;
    count1=1;
    pot_ind=current_ind+N6_ind;
    T_pot_ind=T(pot_ind);
    T_pot_ind(L(pot_ind)~=Labels(i) | KT(pot_ind)~=1)=NaN;
%     if Labels(i)==L(current_ind)
%         [temp_v,temp_ind]=nanmin((T_pot_ind-T(current_ind)).*h1);
%     else
%         [temp_v,temp_ind]=nanmax(-T_pot_ind);
%     end
    [temp_v,temp_ind]=nanmax(-T_pot_ind.*h1);
    while ~isempty(current_ind) && T(current_ind)~=0 && temp_v<=0
        current_ind=pot_ind(temp_ind);
        count1=count1+1;
        count=count+1;
        GDpath(count)=current_ind;
        
        pot_ind=current_ind+N6_ind;
        T_pot_ind=T(pot_ind);
        T_pot_ind(L(pot_ind)~=Labels(i) | KT(pot_ind)~=1)=NaN;
        [temp_v,temp_ind]=nanmin((T_pot_ind-T(current_ind)).*h1);
    end
    
    if T(current_ind)==0
        current_ind=GDpath(count);
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
    else
        disp('Gradient descent failed to find a T=0 voxel.')
        count=count-1;
        End=lengthAM+count; 
        AMlbl(End,End)=0;
        GDpath(count-count1+2:count+1)=NaN;
    end
end

GDpath=GDpath(1:count);

sizeIm=size(Im);
[GDx,GDy,GDz]=ind2sub(sizeIm,GDpath);
min_x=min(GDx); max_x=max(GDx);
min_y=min(GDy); max_y=max(GDy);
min_z=min(GDz); max_z=max(GDz);

pad1=10;
Xmin=max(min_x-1-pad1,0)+1;
Xmax=min(max_x-sizeIm(1)+pad1,0)+sizeIm(1);
Ymin=max(min_y-1-pad1,0)+1;
Ymax=min(max_y-sizeIm(2)+pad1,0)+sizeIm(2);
Zmin=max(min_z-1-pad1,0)+1;
Zmax=min(max_z-sizeIm(3)+pad1,0)+sizeIm(3);

Im_GDpath=Im(Xmin:Xmax,Ymin:Ymax,Zmin:Zmax);
r_GDpath=[GDx-Xmin,GDy-Ymin,GDz-Zmin]+1;
[R_GDpath I_GDpath] = Find_R_fast(Im_GDpath,r_GDpath,R_min,max(0.5,(R_max-R_min)/10),R_max);
I_GDpath(I_GDpath<0)=0;

r_ind=[r_ind;GDpath];
R=[R;R_GDpath];
I=[I,I_GDpath];
GDpath(isnan(GDpath))=[];


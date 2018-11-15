% This function solves of the Eikonal equation by using the Fast Marching algorithm
% of Sethian. T and D are the arival time and distance maps.
% N_steps is the maximum number of re-initialization steps
% Max_dist is the distance at which re-initialization is performed
% Min_dist is the distance to which map is truncated after re-initialization
% SVr contains positions of the seeds

function [AMlbl,r]=FastMarching6(Orig,N_steps,Max_dist,SVr)

background=0.01;
create_loops=true;
pad=2;
Min_dist=inf;
av_window=10;
contrast=0.2;
std_thr=50;

disp('Fast Marching started.')
disp(['# steps  ', '# trees'])

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
clear Orig
sizeIm=sizeOrig+2*pad;
SVr=SVr+pad;
StartVoxel=sub2ind_AS(sizeIm,SVr(:,1),SVr(:,2),SVr(:,3));

AMlbl=sparse(length(StartVoxel),length(StartVoxel));
r_ind=StartVoxel;
[R I] = Find_R_fast(Im,SVr,1,0.5,10);

N6_ind=[-1;+1;-sizeIm(1);+sizeIm(1);-sizeIm(1)*sizeIm(2);+sizeIm(1)*sizeIm(2)];

T=inf(sizeIm); T(StartVoxel)=0;
D=inf(sizeIm); D(StartVoxel)=0;

L=zeros(sizeIm,'uint16'); L(StartVoxel)=(1:length(StartVoxel));
AML=sparse(diag(ones(1,2*length(StartVoxel))));
KT=zeros(sizeIm,'uint8'); KT(StartVoxel)=1; % Known=1, Trial=2
KnownPoints=StartVoxel;

figure(100)
imshow(max(Im,[],3),[0 max(Im(:))]), hold on
[NVx,NVy,NVz]=ind2sub(sizeIm,StartVoxel);
plot3(NVy,NVx,NVz,'r.','MarkerSize',20)
drawnow

stop_condition=false;
Max_T=nan(1,N_steps);
%Max_L=zeros(1,N_steps);
Max_N=nan(1,N_steps);
Max_I=nan(1,N_steps);
M=nan(1,N_steps);
count=0;
Lcount=length(StartVoxel);
pad1=ceil(Max_dist)+pad;
while stop_condition==false && count<=N_steps
    
    % do fast marching untill 1 of 4 conditions is encountered
    [T,D,L,AML,KT,KnownPoints,NewKnownPoint,exit_flag]=Fast_Marching(Im,T,D,L,AML,KT,KnownPoints,N6_ind,Max_dist,inf);
    
    if strcmp(exit_flag,'dist') || strcmp(exit_flag,'collision')
        % find the gradient descent path and relabel the regions
        [T,D,L,AML,KT,AMlbl,r_ind,R,I,GDpath,R_GDpath,I_GDpath,Nregions]=GradientDescent(Im,T,D,L,AML,KT,AMlbl,r_ind,R,I,NewKnownPoint,N6_ind,Max_dist,create_loops);

%         br_ind=find(AMlbl(:,(r_ind==GDpath(end))));
%         R_br=ceil(R(br_ind(1)));
%         C=mean(Im(GDpath(1:end-R_br)))/mean(Im([GDpath(end-R_br+1:end);br_ind(1)]));
        C=mean(I_GDpath(1:ceil(Max_dist/2)))/mean(I);
        
        count=count+1;
        Lcount=Lcount+1-Nregions;
        
        N_trees=2*length(StartVoxel)-(sum(AML(:))-length(AML))/2;
        disp([count,N_trees]);
        
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
            L1=zeros(sizeIm1,'uint16'); L1(GDind1)=L(GDpath(1));
            AML1=sparse(diag(ones(1,2*length(StartVoxel))));
            KT1=zeros(sizeIm1,'uint8'); KT1(GDind1)=1;
            
            Max_T(count)=max(T(KT==1));
            %Max_L(count)=nnz(T==0);
            Max_N(count)=nnz(KT==1);
            Max_I(count)=sum(Im(KT==1));
            [T1,D1,L1,~,~,~,~,~]=Fast_Marching(Im1,T1,D1,L1,AML1,KT1,GDind1,N6_ind1,inf,Max_T(count));
            
            %             Rtemp=sum(D1(KT1==1).*Im1(KT1==1))./sum(Im1(KT1==1));
            %             for jj=1:length(GDpath)
            %                 R(find(r_ind==GDpath(jj)))=Rtemp;
            %             end
            %             br_ind=find(AMlbl(:,find(r_ind==GDpath(end))));
            %             R_br=ceil(R(br_ind(1)));
            %[~,p1]=ttest2(Im1(KT1==1),contrast*Im1(:),0.05,'right');
            
            % merge the maps
            T1=T1(pad+1:end-pad,pad+1:end-pad,pad+1:end-pad);
            D1=D1(pad+1:end-pad,pad+1:end-pad,pad+1:end-pad);
            L1=L1(pad+1:end-pad,pad+1:end-pad,pad+1:end-pad);
            %KT1=KT1(pad+1:end-pad,pad+1:end-pad,pad+1:end-pad);
            
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
            %KT1(~temp_ind & KT(Xmin:Xmax,Ymin:Ymax,Zmin:Zmax)==1)=inf;
            
            temp_ind=(T1<=T(Xmin:Xmax,Ymin:Ymax,Zmin:Zmax));
            tempT=T(Xmin:Xmax,Ymin:Ymax,Zmin:Zmax); tempT(temp_ind)=T1(temp_ind); T(Xmin:Xmax,Ymin:Ymax,Zmin:Zmax)=tempT;
            tempD=D(Xmin:Xmax,Ymin:Ymax,Zmin:Zmax); tempD(temp_ind)=D1(temp_ind); D(Xmin:Xmax,Ymin:Ymax,Zmin:Zmax)=tempD;
            tempL=L(Xmin:Xmax,Ymin:Ymax,Zmin:Zmax); tempL(temp_ind)=L1(temp_ind); L(Xmin:Xmax,Ymin:Ymax,Zmin:Zmax)=tempL;
            %tempKT=KT(Xmin:Xmax,Ymin:Ymax,Zmin:Zmax); tempKT(temp_ind)=KT1(temp_ind); KT(Xmin:Xmax,Ymin:Ymax,Zmin:Zmax)=tempKT;
        end
        
        figure(100)
        [NVx,NVy,NVz]=ind2sub(sizeIm,GDpath);
        plot3(NVy,NVx,NVz,'g.','MarkerSize',6)
        drawnow
        
%         figure(20)
%         subplot(1,2,1)
%         plot(R_GDpath)
%         ylim([0 10])
%         subplot(1,2,2)
%         plot(I_GDpath)
%         ylim([0 1])
         
        delT=Max_T(count)-Max_T(max(1,count-1));
        %delL=Max_L(count)-Max_L(max(1,count-1));
        delN=Max_N(count)-Max_N(max(1,count-1));
        delI=Max_I(count)-Max_I(max(1,count-1));
        
        M(count)=delT*delN*delI;
        MM=(M(count)-nanmean(M(max(1,count-av_window):max(1,count-1))))/nanstd(M(max(1,count-av_window):max(1,count-1)));
        %         figure(10), hold on
        %         plot(count,M(count),'c*')
        disp([MM,C])
        if (MM>std_thr && count>=av_window) || C<contrast
            reply = input('Do you want to continue? ', 's');
            if strcmp(reply,'n')
                stop_condition=true;
            elseif strcmp(reply,'y')
                if (MM>std_thr && count>=av_window)
                    %std_thr=std_thr*1.2;
                    disp(MM)
                end
                if C<contrast
                    %contrast=contrast*0.9;
                    disp(C)
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
    elseif strcmp(exit_flag,'end')
        stop_condition=true;
    end
end

[rx,ry,rz]=ind2sub(sizeIm,r_ind);
r=[rx,ry,rz]-pad;

rem=(sum(AMlbl,1)==0);
AMlbl(rem,:)=[];
AMlbl(:,rem)=[];
r(rem,:)=[];
R(rem,:)=[];

if create_loops
    AMlbl = LabelBranchesAM(AMlbl);
    [AMlbl, r, ~, ~]=Reduce_Small_Loops(AMlbl,r,Max_dist);
end
AMlbl = LabelTreesAM(AMlbl);

disp('Fast Marching is complete.')

% figure(300)
% imshow(max(Im,[],3),[0 max(Im(:))])
% hold on
% PlotAM(AMlbl, r+pad)
% plot3(SVr(:,2),SVr(:,1),SVr(:,3),'r.')
%
% if length(AMlbl)~=size(r,1) || nnz(diag(AMlbl))>0 || nnz(AMlbl-AMlbl')>0
%     disp('Error in AMlbl')
% end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% this function performs Fast Marching untill a stop condition is met
function [T,D,L,AML,KT,KnownPoints,exit_ind,exit_flag]=Fast_Marching(Im,T,D,L,AML,KT,KnownPoints,N6_ind,Max_dist,Max_time)

NHood=[KnownPoints+N6_ind(1);KnownPoints+N6_ind(2);KnownPoints+N6_ind(3); ...
    KnownPoints+N6_ind(4);KnownPoints+N6_ind(5);KnownPoints+N6_ind(6)];
NHood(KT(NHood)==1)=[];
%NHood(KT(NHood)>=1)=[];
NHood=unique(NHood);

TrialPoints=NHood(~isinf(T(NHood))); %!!!!!!!!!!!
T_TrialPoints=T(TrialPoints);

stop_cond=true;
if isempty(NHood)
    stop_cond=false;
    exit_flag='end';
end

while stop_cond
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
    N=ones(size(Tpm_xyz,1),1)*(1:size(Tpm_xyz,2));
    
    ct=1./Im(NHood).^2;
    indt=((cumsum(Tpm_xyz.^2,2)-cumsum(Tpm_xyz,2).^2./N)<=ct*ones(1,size(Tpm_xyz,2)));
    Tpm_xyz(~indt)=0;
    nt=sum(indt,2);
    Tmean=sum(Tpm_xyz,2)./nt;
    T(NHood)=Tmean+(ct./nt-(sum(Tpm_xyz.^2,2)./nt-Tmean.^2)).^0.5;
    
    cd=ones(size(NHood));
    indd=((cumsum(Dpm_xyz.^2,2)-cumsum(Dpm_xyz,2).^2./N)<=cd*ones(1,size(Dpm_xyz,2)));
    Dpm_xyz(~indd)=0;
    nd=sum(indd,2);
    Dmean=sum(Dpm_xyz,2)./nd;
    D(NHood)=Dmean+(cd./nd-(sum(Dpm_xyz.^2,2)./nd-Dmean.^2)).^0.5;
    
    temp=(1:length(NHood))'+length(NHood).*(IX(:,1)-1);
    L(NHood)=L(ind(temp));
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    keep=NHood(KT(NHood)==0 & ~isinf(T(NHood))); %!!!!!!!!!!!
    KT(keep)=2;
    TrialPoints=[TrialPoints;keep];
    T_TrialPoints=[T_TrialPoints;T(keep)];
    
    [~,min_ind]=min(T_TrialPoints);
    NewKnownPoint=TrialPoints(min_ind);
    
    %     if length(TrialPoints)<100
    %         1
    %     end
    if ~isinf(T(NewKnownPoint)) && ~isempty(NewKnownPoint)
        KT(NewKnownPoint)=1;
        T_TrialPoints(min_ind)=inf;
        NHood=NewKnownPoint+N6_ind;
        KnownPoints=[KnownPoints;NewKnownPoint];
        
        tempL=L(NHood(KT(NHood)==1));
        collision=(sum(AML(L(NewKnownPoint),tempL))~=length(tempL));
        
        if collision==1
            AML(L(NewKnownPoint),tempL)=1;
            AML(tempL,L(NewKnownPoint))=1;
        end
        NHood(KT(NHood)==1)=[];
        %NHood(KT(NHood)>=1)=[];
        
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
% This function performs gradient descent from ind (a single point) to
% T=0. T map of the region is reset to the gradient descent path, Tnew.
% If there are multiple labels neighboring ind, the gradient descent will
% go into all these labeld regions. The result will be merged and
% relabeled.
function [T,D,L,AML,KT,AMlbl,r_ind,R,I,GDpath,R_GDpath,I_GDpath,Nregions]=GradientDescent(Im,T,D,L,AML,KT,AMlbl,r_ind,R,I,ind,N6_ind,Max_dist,create_loops)

lengthAM=size(AMlbl,1);
pot_ind=ind+N6_ind;
pot_ind(KT(pot_ind)~=1)=[];
Labels=double(unique(L(pot_ind)));
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
[R_GDpath I_GDpath] = Find_R_fast(Im_GDpath,r_GDpath,1,0.5,10);

r_ind=[r_ind;GDpath];
R=[R;R_GDpath];
I=[I,I_GDpath];
GDpath(isnan(GDpath))=[];

T(GDpath)=0;
D(GDpath)=0;
KT(GDpath)=1;

% relabel
if create_loops
    L_new=double(max(L(:))+1);
    AML(L_new,Labels)=1;
    AML(Labels,L_new)=1;
elseif ~create_loops
    if Nregions>1
        L_new=max(L(:))+1;
        for i=1:Nregions
            L((L==Labels(i)))=L_new;
            AMlbl((AMlbl==Labels(i)))=L_new;
        end
    else
        L_new=L(ind);
    end
    L(GDpath)=L_new;
end
AMlbl(isnan(AMlbl))=L_new;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%This function locates tips and extends them in a direction determined by
%coneangle, by a maximum length given by maxextension.
function [AMlbl,r] = ExtendTips(Orig,AMlbl,r,Max_dist,unisotropy)

unisotropy(unisotropy<0)=0;
if sum(unisotropy)==0
    unisotropy=[1,1,1];
end
unisotropy=unisotropy./sum(unisotropy.^2)^0.5;
h1=[unisotropy(1),unisotropy(1),unisotropy(2),unisotropy(2),unisotropy(3),unisotropy(3)]';
h2=[unisotropy(1),unisotropy(1),unisotropy(2),unisotropy(2),unisotropy(3),unisotropy(3)].^2;

%Parameter list
updatevol=Max_dist;
maxextension=ceil(Max_dist/2);
coneangle=1/3^0.5; %[0,0,0] to [1,1,1] in 3D
background=0.01;

Orig=double(Orig);
Orig=Orig-min(Orig(:));
Orig=(Orig./max(Orig(:)));
Orig(Orig<background)=background;
pad=2;
Orig=padarray(Orig,[pad,pad,pad],0);
r=r+pad;
sizeIm=size(Orig);

r_ind=sub2ind(sizeIm,r(:,1),r(:,2),r(:,3));
Tmap=inf(size(Orig));
Tmap(r_ind)=0;

%Locate tips
tips=StepBackTips(AMlbl,Max_dist);
tipstartind=zeros(length(tips),1);
tipendind=tipstartind;
L_tips=tipstartind;
for i=1:size(tips,1)
    tipstartind(i)=tips{i}(1);
    tipendind(i)=tips{i}(end);
    L_tips(i)=AMlbl(tips{i}(1),tips{i}(2));
end
tipstartind=r_ind(tipstartind);
tipendind=r_ind(tipendind);
tipstartsub=zeros(length(tipstartind),3);
tipendsub=zeros(length(tipstartind),3);
[tipstartsub(:,1),tipstartsub(:,2),tipstartsub(:,3)]=ind2sub(sizeIm,tipstartind);
[tipendsub(:,1),tipendsub(:,2),tipendsub(:,3)]=ind2sub(sizeIm,tipendind);

%Find normalized tip direction vector
nvec=tipstartsub-tipendsub;
nvecnorm=(nvec(:,1).^2+nvec(:,2).^2+nvec(:,3).^2).^(0.5);
nvec=nvec./repmat(nvecnorm,1,3);

% Enter for loop for individual tip
for i=1:length(tipstartind)
    
    %Volume selection
    Xmin=max(tipstartsub(i,1)-1-updatevol,0)+pad+1;
    Xmax=min(tipstartsub(i,1)-sizeIm(1)+updatevol,0)+sizeIm(1)-pad;
    Ymin=max(tipstartsub(i,2)-1-updatevol,0)+pad+1;
    Ymax=min(tipstartsub(i,2)-sizeIm(2)+updatevol,0)+sizeIm(2)-pad;
    Zmin=max(tipstartsub(i,3)-1-updatevol,0)+pad+1;
    Zmax=min(tipstartsub(i,3)-sizeIm(3)+updatevol,0)+sizeIm(3)-pad;
    
    %Initializations
    Im1=Orig(Xmin:Xmax,Ymin:Ymax,Zmin:Zmax);
    Im1=padarray(Im1,[pad,pad,pad],0);
    T1=Tmap(Xmin:Xmax,Ymin:Ymax,Zmin:Zmax);
    T1=padarray(T1,[pad,pad,pad],inf);
    sizeIm1=size(Im1);
    N6_ind1=[-1;+1;-sizeIm1(1);+sizeIm1(1);-sizeIm1(1)*sizeIm1(2);+sizeIm1(1)*sizeIm1(2)];
    
    %Expres indices of points relative to local volume
    Tipstartx1=tipstartsub(i,1)-Xmin+1+pad;
    Tipstarty1=tipstartsub(i,2)-Ymin+1+pad;
    Tipstartz1=tipstartsub(i,3)-Zmin+1+pad;
    Tipstartind1=sub2ind(sizeIm1,Tipstartx1,Tipstarty1,Tipstartz1);
    
    %Initialize other matrices for fastMarching on local volume
    T1(T1>0)=inf; T1(Tipstartind1)=0;
    KnownPoints1=find(T1==0);
    
    %Keep KnownPoints1 only if they belong to the cone
    [NHx,NHy,NHz]=ind2sub(sizeIm1,KnownPoints1);
    NHrx=NHx-Tipstartx1;
    NHry=NHy-Tipstarty1;
    NHrz=NHz-Tipstartz1;
    KeepTrialind1=(NHrx==0 & NHry==0 & NHrz==0);
    KeepTrialind2=(nvec(i,1)*NHrx+nvec(i,2)*NHry+nvec(i,3)*NHrz)./(NHrx.^2+NHry.^2+NHrz.^2).^.5;
    KnownPoints1=KnownPoints1(KeepTrialind1 | KeepTrialind2>=coneangle ); 
    
    %Use KnownPoints1 to initialize other matrices
    D1=inf(sizeIm1); D1(KnownPoints1)=0;
    L1=zeros(sizeIm1); L1(T1==0)=L_tips(i)+1; L1(Tipstartind1)=L_tips(i);
    KT1=zeros(sizeIm1); KT1(KnownPoints1)=1;
    
    %Modified FastMarching code begins
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %NHood1 may contain points belonging to pad.
    NHood1=[KnownPoints1+N6_ind1(1);KnownPoints1+N6_ind1(2);KnownPoints1+N6_ind1(3); ...
        KnownPoints1+N6_ind1(4);KnownPoints1+N6_ind1(5);KnownPoints1+N6_ind1(6)];
    NHood1(KT1(NHood1)==1)=[];
    NHood1=unique(NHood1);
    
    TrialPoints1=[];
    T_TrialPoints1=[];
    
    stop_cond=true;
    while stop_cond
        if ~isempty(NHood1)
            ind=[NHood1+N6_ind1(1),NHood1+N6_ind1(2),NHood1+N6_ind1(3),NHood1+N6_ind1(4),NHood1+N6_ind1(5),NHood1+N6_ind1(6)];
            
            %Backsolver begins
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            KT_1=(KT1(ind)~=1);
            
            Tpm_xyz=T1(ind);
            Tpm_xyz(KT_1)=inf;
            [Tpm_xyz,IX]=sort(Tpm_xyz,2);
            
            Dpm_xyz=D1(ind);
            Dpm_xyz(KT_1)=inf;
            Dpm_xyz=sort(Dpm_xyz,2);
            
            H=h2(IX);
            H_cum=cumsum(H,2);
            ct=1./Im1(NHood1).^2;
            Tpm_xyz_cum=cumsum(Tpm_xyz.*H,2);
            Tpm_xyz2_cum=cumsum(Tpm_xyz.^2.*H,2);
            nt=sum(((Tpm_xyz2_cum-Tpm_xyz_cum.^2./H_cum)<=ct*ones(1,size(Tpm_xyz,2))),2);
            if sum(h2==0)>0
                nt(nt==0)=1;
            end
            ind_nt=(1:size(Tpm_xyz,1))'+(nt-1).*size(Tpm_xyz,1);
            ntH=H_cum(ind_nt);
            temp=Tpm_xyz_cum(ind_nt);
            T1(NHood1)=(temp+(ct.*ntH-(Tpm_xyz2_cum(ind_nt).*ntH-temp.^2)).^0.5)./ntH;
            
            N=ones(size(Tpm_xyz,1),1)*(1:size(Tpm_xyz,2));
            Dpm_xyz_cum=cumsum(Dpm_xyz,2);
            Dpm_xyz2_cum=cumsum(Dpm_xyz.^2,2);
            nd=sum(((Dpm_xyz2_cum-Dpm_xyz_cum.^2./N)<=ones(size(Dpm_xyz))),2);
            ind_nd=(1:size(Dpm_xyz,1))'+(nd-1).*size(Dpm_xyz,1);
            temp=Dpm_xyz_cum(ind_nd);
            D1(NHood1)=(temp+(nd-(Dpm_xyz2_cum(ind_nd).*nd-temp.^2)).^0.5)./nd;
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            temp=(1:length(NHood1))'+length(NHood1).*(IX(:,1)-1);
            L1(NHood1)=L1(ind(temp));
            
            keep=NHood1(KT1(NHood1)==0);
            KT1(keep)=2;
            TrialPoints1=[TrialPoints1;keep];
            T_TrialPoints1=[T_TrialPoints1;T1(keep)];
            
        end
        
        [T_min,min_ind]=min(T_TrialPoints1);
        NewKnownPoint1=TrialPoints1(min_ind);
        if ~isempty(NewKnownPoint1) && T_min<inf
            KnownPoints1=[KnownPoints1;NewKnownPoint1];
            KT1(NewKnownPoint1)=1;
            T_TrialPoints1(min_ind)=inf;
            NHood1=NewKnownPoint1+N6_ind1;
           
            tempL=L1(NHood1(KT1(NHood1)==1));
            tempL=tempL(tempL>0);
            collision=(sum(tempL-tempL(1))~=0);
            NHood1(KT1(NHood1)==1)=[];
            
            %Exit conditions
            if D1(NewKnownPoint1)>=ceil(maxextension) && L1(NewKnownPoint1)==L_tips(i)
                stop_cond=false;
            end
            
            if collision
                stop_cond=false;
                T1(L1~=L_tips(i))=inf;
            end
            
            %Refresh NHood for next iteration: points in the backward direction
            NHood1=unique(NHood1);
            [NHx,NHy,NHz]=ind2sub(sizeIm1,NHood1);
            NHrx=NHx-Tipstartx1;
            NHry=NHy-Tipstarty1;
            NHrz=NHz-Tipstartz1;
            KeepTrialind1=(NHrx==0 & NHry==0 & NHrz==0);
            KeepTrialind2=(nvec(i,1)*NHrx+nvec(i,2)*NHry+nvec(i,3)*NHrz)./(NHrx.^2+NHry.^2+NHrz.^2).^.5;
            NHood1=NHood1(KeepTrialind1 | KeepTrialind2>=coneangle );
            
        else
            NewKnownPoint1=[];
            stop_cond=false;
        end
    end
    exit_ind=NewKnownPoint1;
    
    %Modified Gradient Descent begins
    if ~isempty(exit_ind)
        lengthAM=size(AMlbl,1);
        GDpath1=zeros(2*maxextension,1);
        GDpath1(1)=exit_ind;
        AMlbl(lengthAM+1,lengthAM+1)=0;
        count=1;
        
        current_ind=exit_ind;
        prev_ind=current_ind;
        T_temp=0.1;
        while T1(current_ind)~=0 && ~isempty(current_ind) && T_temp>0
            pot_ind1=current_ind+N6_ind1;
            temp=L1(pot_ind1)~=L_tips(i) | KT1(pot_ind1)~=1 | pot_ind1==prev_ind;
            pot_ind1(temp)=[];
            h1_temp=h1(~temp);
            
            [T_temp,temp]=nanmax((T1(current_ind)-T1(pot_ind1)).*h1_temp);
            prev_ind=current_ind;
            current_ind=pot_ind1(temp);
            
            if ~isempty(current_ind) 
                
                count=count+1;
                GDpath1(count)=current_ind;
            end
        end
        
        GDpath1=GDpath1(1:count);
        [GDx1,GDy1,GDz1]=ind2sub(sizeIm1,GDpath1);
        GDx=GDx1+Xmin-1-pad;
        GDy=GDy1+Ymin-1-pad;
        GDz=GDz1+Zmin-1-pad;
        GDpath=sub2ind(sizeIm,GDx,GDy,GDz);
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
        
        GDpath=GDpath(1:count);
        r_ind=[r_ind;GDpath];
        AMlbl(isnan(AMlbl))=L_tips(i);
    end
end
[rx,ry,rz]=ind2sub(sizeIm,r_ind);
r=[rx,ry,rz];
r=r-pad;



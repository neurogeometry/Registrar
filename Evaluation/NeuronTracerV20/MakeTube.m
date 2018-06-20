function [] = MakeTube(AM,r,nsteps,zmod)
%This function generates a binary tube in 3d around the input trace.

cyl_size=50; %Size of ~cylinder drawn around the trace
Ithr=1000; %Change threshold after figures are plotted using 'caxis([0 new_threshold_value])'
%This script creates a binary filter around the linear axon traces and uses
%it to obtain z-projection of the filtered image.
animal = {'DL088'};a=1;
section = {'001'};s=1;
timepoint = {'I'};
axon = [6];
channel = {'G'};

im_pth=ismacispc('C:\Users\Rohan\Desktop\DatasetIM\');
manualtrace_pth=ismacispc('C:\Users\Rohan\Dropbox\Lab\Plasticity\dat\OptimizedTraces\');
%Expected file on manualtrace_pth:
%C:\Users\Rohan\Dropbox\Lab\Plasticity\dat\ManualTraces\DL088B005\A001.swc

for t=1:length(timepoint)
    IM=load(([im_pth,animal{a},timepoint{t},section{s},channel{1}]),'Original');
    IM=IM.Original;
    sizeIm=size(IM);
    AM=[];r=[];I=[];NHood=[];
    for ax=1:length(axon)
        axonstr=sprintf('A%.3f',10^-3*axon(ax));
        axonstr(2:3)=[];
        loadfile=(ismacispc([manualtrace_pth,animal{a},timepoint{t},section{s},'\',axonstr]));
        [AM,r]=SWC2AM(([loadfile,'.swc']));
        R=zeros(size(r,1),1);
        [~,rceil,~] = AdjustPPM(AM,r,R,1);
        rceil=ceil(rceil);
        rceil((rceil(:,1)>sizeIm(1) | rceil(:,2)>sizeIm(2) | rceil(:,3)>sizeIm(3)),:)=[];
        NHood=[NHood;unique(sub2ind(sizeIm,rceil(:,1),rceil(:,2),rceil(:,3)))];
    end
    
    N6_ind=[-1;+1;-sizeIm(1);+sizeIm(1);-sizeIm(1)*sizeIm(2);+sizeIm(1)*sizeIm(2)];
    for i=1:cyl_size
        NHood=[NHood+N6_ind(1);NHood+N6_ind(2);NHood+N6_ind(3); ...
            NHood+N6_ind(4)];
        NHood(NHood<1 | NHood>numel(IM))=[];
        NHood=unique(NHood);
        if mod(i,zmod)==0
            NHood=[NHood+N6_ind(1);NHood+N6_ind(2);NHood+N6_ind(3); ...
                NHood+N6_ind(4);NHood+N6_ind(5);NHood+N6_ind(6)];
            NHood(NHood<1 | NHood>numel(IM))=[];
            NHood=unique(NHood);
        end
    end
    
    IMfig=zeros(size(IM),class(IM));
    IMfig(NHood)=IM(NHood);
    
    minr=round(min(r,[],1));
    maxr=round(max(r,[],1));
    
    minx=max(minr(1)-10,1);
    miny=max(minr(2)-10,1);
    minz=max(minr(3)-10,1);
    maxx=min(maxr(1)+10,size(IMfig,1));
    maxy=min(maxr(2)+10,size(IMfig,2));
    maxz=min(maxr(3)+10,size(IMfig,3));
    
    figure(2),clf(2),imshow(squeeze(max(IMfig(minx:maxx,miny:maxy,minz:maxz),[],3)),[0 Ithr]),hold on
    PlotAM_c(AM,[r(:,1)-minx+1,r(:,2)-miny+1,r(:,3)-minz+1],[1 0 0]);
    title('Z Projection')
    % movegui(1,'west')
    %
    %     figure(2),clf(2),imshow(squeeze(max(IMfig(minx:maxx,miny:maxy,minz:maxz),[],2)),[0 Ithr]),hold on
    %     PlotAM_c(AM,[r(:,1)-minx,r(:,3)-minz,r(:,2)-miny],[1 0 0])
    %     title('Y Projection')
    %     movegui(2,'center')
    %
    %     figure(3),clf(3),imshow(squeeze(max(IMfig(minx:maxx,miny:maxy,minz:maxz),[],1)),[0 Ithr]),hold on
    %     PlotAM_c(AM,[r(:,2)-miny,r(:,3)-minz,r(:,1)-minx],[1 0 0]);
    %     title('X Projection')
    %     movegui(3,'east')
    
    drawnow
    display(['Time: ',timepoint{t},  ' Axon: ',axonstr])
end

%%%%%%%%%%%%%%%%%%%%%%%%%%






end

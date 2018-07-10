function showNRRegResult (StackList,sourceID,targetID,Boutons,Traces,showpairs)
N=1;
Source_Stack_File = char(StackList(sourceID,1));
Target_Stack_File = char(StackList(targetID,1));
IM_Source=ImportStack(char(Source_Stack_File));
IM_Source = uint8(double(IM_Source)./double(max(IM_Source(:))).*255);
IM_Target=ImportStack(char(Target_Stack_File));
IM_Target = uint8(double(IM_Target)./double(max(IM_Target(:))).*255);
IM_source_max=max(IM_Source,[],3);
IM_target_max=max(IM_Target,[],3);

if showpairs
%% Before Registration
N=N+1;
figure(N),imshowpair(IM_source_max,IM_target_max,'Scaling','independent')
title('Before Registration');


%% After Registration
[IM_Target_NR,StackPosition_prime,~]=Perform_Bspline_Transform(IM_Target,[1;1;1],L,b,Cxyz,Nxyz,nxyz,Grid_start,affine);
IM_Target_NR_max=max(IM_Target_NR,[],3);
MIN=min([1;1],StackPosition_prime(1:2));
MAX=max(size(IM_source_max)',size(IM_target_max)'+StackPosition_prime(1:2)-1);
temp=zeros(MAX(1)-MIN(1)+1,MAX(2)-MIN(2)+1,'uint8');

IM_Target_NR_max_P=temp;
IM_Target_NR_max_P(StackPosition_prime(1)-MIN(1)+1:StackPosition_prime(1)-MIN(1)+size(IM_Target_NR_max,1),...
    StackPosition_prime(2)-MIN(2)+1:StackPosition_prime(2)-MIN(2)+size(IM_Target_NR_max,2))=IM_Target_NR_max;
IM_source_max_P=temp;
IM_source_max_P(2-MIN(1):1-MIN(1)+size(IM_source_max,1),...
    2-MIN(2):1-MIN(2)+size(IM_source_max,2))=IM_source_max;
N=N+1;
figure(N),imshowpair(IM_Target_NR_max_P,IM_source_max_P,'Scaling','independent')
end

%% Boutons
if ~isempty(Boutons) 
    N=N+1;
    figure(N),imshow(IM_target_max,[0 20]);
    hold on; plot(Boutons.r2(:,2),Boutons.r2(:,1),'*r')
    hold on; plot(Boutons.r1(:,2),Boutons.r1(:,1),'*g')
    title('Boutons Before Registration');
    N=N+1;
    figure(N),imshow(IM_target_max,[0 20]);
    hold on; plot(Boutons.r2(:,2),Boutons.r2(:,1),'*r')
    hold on; plot(Boutons.r1_NR(:,2),Boutons.r1_NR(:,1),'*g')
    title('Boutons After Registration');
%     N=N+1;
%     figure(N),imshow(IM_source_max,[0 20]);
%     hold on; plot(Boutons.r1(:,2),Boutons.r1(:,1),'*r')
%     hold on; plot(Boutons.r1_NR(:,2),Boutons.r1_NR(:,1),'*g')
%     title('Boutons Movement After Registration');
end
%% Traces
if ~isempty(Traces)
    figure;imshow(IM_source_max,[0 500]);
    hold on; PlotAM(AM_Source,r_Source,'r')
    figure;imshow(IM_target_max,[0 500]);
    hold on; PlotAM(AM_Target,r_Target,'g')

end









































%%

%                     IM_source_maxR = imresize(IM_source_max,size(IM_Target_NR_max));
%                     figure(3),imshowpair(IM_Target_NR_max,IM_source_maxR,'Scaling','independent')

%         [IM_Source_NR,~]=Perform_Bspline_Transform(IM_Source,[0,0,0],L,b,Cxyz,Nxyz,nxyz,Grid_start,affine);
%         IM_Source_NR_max=max(IM_Source_NR,[],3);
%         figure(3)
%         imshowpair(IM_Source_NR_max,IM_target_max,'Scaling','independent')
%
%         [optimizer, metric] = imregconfig('multimodal');
%         optimizer.InitialRadius = 0.009;
%         optimizer.Epsilon = 1.5e-4;
%         optimizer.GrowthFactor = 1.01;
%         optimizer.MaximumIterations = 300;
%         movingRegistered = imregister(IM_Source, IM_Target, 'affine', optimizer, metric);
%         figure(4)
%         imshowpair(IM_target_max,max(movingRegistered,[],3),'Scaling','independent')

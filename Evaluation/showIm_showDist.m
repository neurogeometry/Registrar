filename = '..\..\RegistrationEvaluation\result_NR_fiji';
load([filename,'.mat'])
load([filename,'\Evaluation.mat'])
tic
Functionspath = '..\';
addpath([Functionspath,'Functions']);
showDiff = 0;
csvPath = '..\..\RegistrationEvaluation\TimeLapse_Holtmaat_StackList.csv';
StackList = table2cell(readtable(csvPath,'Delimiter',','));

% IM_source_maxproj = cell(1,size(result,2));
% IM_Target_NR_maxproj = cell(1,size(result,2));

if showDiff
    for sourceID = 1:size(result,2)
        targetID = sourceID + 1;
        
        Source_Stack_File = char(StackList(sourceID,1));
        Target_Stack_File = char(StackList(targetID,1));
        IM_Source=ImportStack(char(Source_Stack_File));
        IM_Source = uint8(double(IM_Source)./double(max(IM_Source(:))).*255);
        IM_Target=ImportStack(char(Target_Stack_File));
        IM_Target = uint8(double(IM_Target)./double(max(IM_Target(:))).*255);
        IM_source_max=max(IM_Source,[],3);
        IM_target_max=max(IM_Target,[],3);
        clear IM_Source
        
        
        
        %before
        %         figure(2*sourceID-1)
        %         imshowpair(IM_source_max,IM_target_max,'Scaling','independent')
        %
        %                     (1:500,1:500,1:100)
        
        [IM_Target_NR,StackPosition_prime,~]=Perform_Bspline_Transform(IM_Target,[1;1;1],result{1,sourceID}.T.L,result{1,sourceID}.T.b,result{1,sourceID}.T.Cxyz,result{1,sourceID}.T.Nxyz,result{1,sourceID}.T.nxyz,result{1,sourceID}.T.Grid_start,result{1,sourceID}.T.affine);
        clear IM_Target
        IM_Target_NR_max=max(IM_Target_NR,[],3);
        clear IM_Target_NR
        MIN=min([1;1],StackPosition_prime(1:2));
        MAX=max(size(IM_source_max)',size(IM_target_max)'+StackPosition_prime(1:2)-1);
        temp=zeros(MAX(1)-MIN(1)+1,MAX(2)-MIN(2)+1,'uint8');
        
        IM_Target_NR_max_P=temp;
        IM_Target_NR_max_P(StackPosition_prime(1)-MIN(1)+1:StackPosition_prime(1)-MIN(1)+size(IM_Target_NR_max,1),...
            StackPosition_prime(2)-MIN(2)+1:StackPosition_prime(2)-MIN(2)+size(IM_Target_NR_max,2))=IM_Target_NR_max;
        IM_source_max_P=temp;
        IM_source_max_P(2-MIN(1):1-MIN(1)+size(IM_source_max,1),...
            2-MIN(2):1-MIN(2)+size(IM_source_max,2))=IM_source_max;
        
        %         figure(2*sourceID),imshowpair(IM_Target_NR_max_P,IM_source_max_P,'Scaling','independent')
        save([filename,'\IMs',num2str(sourceID),'.mat'],'IM_Target_NR_max_P','IM_source_max_P')
        disp(sourceID)
    end
end

hist1 = hist(All_Dis_voxel,[0.5:1:29.5]);
hist2 = hist(All_Dis_NonRigid_voxel,[0.5:1:29.5]);
hist3 = hist(All_Dis_fiji_Affine_voxel,[0.5:1:29.5]);

figure
plot(hist1)
hold on
plot(hist2)
hold on
plot(hist3)
% legend(['trace mismatch before registration:','mean=',num2str(mean(All_Dis_um)),',std=',num2str(std(All_Dis_um))],['trace mismatch after registration:','mean=',num2str(mean(All_Dis_NonRigid_um)),',std=',num2str(std(All_Dis_NonRigid_um))])
legend(['trace mismatch before registration:','mean=',num2str(mean(All_Dis_voxel)),',std=',num2str(std(All_Dis_voxel))],['trace mismatch after registration:','mean=',num2str(mean(All_Dis_NonRigid_voxel)),',std=',num2str(std(All_Dis_NonRigid_voxel))],['trace mismatch after registration:','mean=',num2str(mean(All_Dis_fiji_Affine_voxel)),',std=',num2str(std(All_Dis_fiji_Affine_voxel))])


hist1 = hist(All_Dis_um,[0.5:1:29.5]);
hist2 = hist(All_Dis_NonRigid_um,[0.5:1:29.5]);
hist3 = hist(All_Dis_fiji_Affine_um,[0.5:1:29.5]);


figure
plot(hist1)
hold on
plot(hist2)
hold on
plot(hist3)
% legend(['trace mismatch before registration:','mean=',num2str(mean(All_Dis_um)),',std=',num2str(std(All_Dis_um))],['trace mismatch after registration:','mean=',num2str(mean(All_Dis_NonRigid_um)),',std=',num2str(std(All_Dis_NonRigid_um))])
legend(['trace mismatch before registration:','mean=',num2str(mean(All_Dis_voxel)),',std=',num2str(std(All_Dis_voxel))],['trace mismatch after registration:','mean=',num2str(mean(All_Dis_NonRigid_voxel)),',std=',num2str(std(All_Dis_NonRigid_voxel))],['trace mismatch after registration:','mean=',num2str(mean(All_Dis_fiji_Affine_voxel)),',std=',num2str(std(All_Dis_fiji_Affine_voxel))])

figure
for i = 1:size(result,2)
    d_bouton = result{1,i}.Bouton.r1-result{1,i}.Bouton.r2;
    plot(d_bouton(:,1),d_bouton(:,2),'.','color','r')
    hold on
    d_bouton_NR = result{1,i}.Bouton.r1_NR-result{1,i}.Bouton.r2;
    plot(d_bouton_NR(:,1),d_bouton_NR(:,2),'.','color','g')
    hold on
end
legend('before registration','after registration')

% save([filename,'\Evaluation.mat'],'All_Dis_um','All_Dis_voxel','All_Dis_NonRigid_um','All_Dis_NonRigid_voxel','eucl_bouton','eucl_bouton_NR')
bhist1 = hist(eucl_bouton,[0.5:1:39.5]);
bhist2 = hist(eucl_bouton_NR,[0.5:1:39.5]);
bhist3 = hist(eucl_bouton_fiji_Affine,[0.5:1:39.5]);
figure,plot(bhist1)
hold on
plot(bhist2)
hold on
plot(bhist3)
legend(['trace mismatch before registration:','mean=',num2str(mean(All_Dis_voxel)),',std=',num2str(std(All_Dis_voxel))],['trace mismatch after registration:','mean=',num2str(mean(All_Dis_NonRigid_voxel)),',std=',num2str(std(All_Dis_NonRigid_voxel))],['trace mismatch after registration:','mean=',num2str(mean(All_Dis_fiji_Affine_voxel)),',std=',num2str(std(All_Dis_fiji_Affine_voxel))])

toc
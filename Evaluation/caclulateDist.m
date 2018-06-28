filename = '..\..\RegistrationEvaluation\result';
load([filename,'.mat'])
mkdir(filename)
tic
All_Dis_um = [];
All_Dis_voxel = [];
All_Dis_NonRigid_um = [];
All_Dis_NonRigid_voxel = [];
for i = 1:size(result,2)
    %     for i = 1:1
    for j = 1:size(result{1,i}.Trace.AM1,2)
        disp([i,j])
        [Dis_um,Dis_voxel] = TraceDistance(result{1,i}.Trace.AM1{j}, result{1,i}.Trace.r1{j}, result{1,i}.Trace.AM2{j}, result{1,i}.Trace.r2{j},pixelSize,0);
        [Dis_NonRigid_um,Dis_NonRigid_voxel] = TraceDistance(result{1,i}.Trace.AM1{j}, result{1,i}.Trace.r1_NR{j}, result{1,i}.Trace.AM2{j}, result{1,i}.Trace.r2{j},pixelSize,0);
        All_Dis_um = [All_Dis_um,Dis_um];
        All_Dis_voxel = [All_Dis_voxel,Dis_voxel];
        All_Dis_NonRigid_um = [All_Dis_NonRigid_um,Dis_NonRigid_um];
        All_Dis_NonRigid_voxel = [All_Dis_NonRigid_voxel,Dis_NonRigid_voxel];
    end
end
hist1 = hist(All_Dis_voxel,[0.5:1:29.5]);
hist2 = hist(All_Dis_NonRigid_voxel,[0.5:1:29.5]);
figure,plot(hist1)
hold on
plot(hist2)
legend(['trace mismatch before registration:','mean=',num2str(mean(All_Dis_um)),',std=',num2str(std(All_Dis_um))],['trace mismatch after registration:','mean=',num2str(mean(All_Dis_NonRigid_um)),',std=',num2str(std(All_Dis_NonRigid_um))])
% % j = 3;
% % figure,PlotAM(result{1,1}.Trace.AM{j}, result{1,1}.Trace.r2{j},'b')
% % hold on
% % PlotAM(result{1,1}.Trace.AM{j}, result{1,1}.Trace.r1{j},'r')
% % hold on
% % PlotAM(result{1,1}.Trace.AM{j}, result{1,1}.Trace.r1_NR{j},'g')
% % legend('Target(b)','Source(r)','Non-Rigid(g)')
eucl_bouton = [];
eucl_bouton_NR = [];
figure
for i = 1:size(result,2)
    d_bouton = result{1,i}.Bouton.r1-result{1,i}.Bouton.r2;
    plot(d_bouton(:,1),d_bouton(:,2),'.','color','r')
    hold on
    d_bouton_NR = result{1,i}.Bouton.r1_NR-result{1,i}.Bouton.r2;
    plot(d_bouton_NR(:,1),d_bouton_NR(:,2),'.','color','g')
    hold on
    eucl_bouton = [eucl_bouton;sqrt(sum(d_bouton.^2,2))];
    eucl_bouton_NR = [eucl_bouton_NR;sqrt(sum(d_bouton_NR.^2,2))];
end
legend('before registration','after registration')
save([filename,'\Evaluation.mat'],'All_Dis_um','All_Dis_voxel','All_Dis_NonRigid_um','All_Dis_NonRigid_voxel','eucl_bouton','eucl_bouton_NR')
bhist1 = hist(eucl_bouton,[0.5:1:39.5]);
bhist2 = hist(eucl_bouton_NR,[0.5:1:39.5]);
figure,plot(bhist1)
hold on
plot(bhist2)
legend(['bouton mismatch before registration:','mean=',num2str(mean(eucl_bouton)),',std=',num2str(std(eucl_bouton))],['bouton mismatch after registration:','mean=',num2str(mean(eucl_bouton_NR)),',std=',num2str(std(eucl_bouton_NR))])
toc
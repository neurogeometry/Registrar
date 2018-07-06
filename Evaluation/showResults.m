clear all;
close all;
mu_Names = {'0','0.00097656','0.0078125','0.0039063','0.0039063','0.015625','0.03125','0.0625','0.125','0.25','0.5','1','2','4','8','16','32','64','128','256','512','1024'};

Calculate_Trace = 1;
NR_mean_Trace = zeros(size(mu_Names,2),1);
NR_std_Trace = zeros(size(mu_Names,2),1);

for nummu = 1:size(mu_Names,2)
    filename = ['..\..\RegistrationEvaluation\MatchedPoints_Non-Rigid_Seyed\result_NR_fiji_',mu_Names{nummu}];
    load([filename,'\Evaluation.mat'])
    NR_mean_Trace(nummu) = mean(All_Dis_NonRigid_voxel);
    NR_std_Trace(nummu) = std(All_Dis_NonRigid_voxel);

%             hist1 = hist(All_Dis_voxel,[0.5:1:29.5]);
%         hist2 = hist(All_Dis_NonRigid_voxel,[0.5:1:29.5]);
%         hist3 = hist(All_Dis_fiji_Affine_voxel,[0.5:1:29.5]);
%         figure,plot(hist1)
%         hold on
%         plot(hist2)
%         hold on
%         plot(hist3)
%         legend(['trace mismatch before registration:','mean=',num2str(mean(All_Dis_voxel)),',std=',num2str(std(All_Dis_voxel))],['trace mismatch after Non-Rigid registration:','mean=',num2str(mean(All_Dis_NonRigid_voxel)),',std=',num2str(std(All_Dis_NonRigid_voxel))],['trace mismatch after fiji-Affine registration:','mean=',num2str(mean(All_Dis_fiji_Affine_voxel)),',std=',num2str(std(All_Dis_fiji_Affine_voxel))])
%
    
%     tic
%     if Calculate_Trace == 1
%         All_Dis_um = [];
%         All_Dis_voxel = [];
%         All_Dis_NonRigid_um = [];
%         All_Dis_NonRigid_voxel = [];
%         All_Dis_fiji_Affine_um = [];
%         All_Dis_fiji_Affine_voxel = [];
%         for i = 1:size(result,2)
% %                     for i = 1:2
% %             for j = 1:size(result{1,i}.Trace.AM1,2)
%                         for j = 1:6
%                 disp([i,j])
%                 [Dis_um,Dis_voxel] = TraceDistance(result{1,i}.Trace.AM1{j}, result{1,i}.Trace.r1{j}, result{1,i}.Trace.AM2{j}, result{1,i}.Trace.r2{j},pixelSize,0);
%                 [Dis_NonRigid_um,Dis_NonRigid_voxel] = TraceDistance(result{1,i}.Trace.AM1{j}, result{1,i}.Trace.r1_NR{j}, result{1,i}.Trace.AM2{j}, result{1,i}.Trace.r2{j},pixelSize,0);
%                 [Dis_fiji_Affine_um,Dis_fiji_Affine_voxel] = TraceDistance(result{1,i}.Trace.AM1{j}, result{1,i}.Trace.r1_fiji_Affine{j}, result{1,i}.Trace.AM2{j}, result{1,i}.Trace.r2_fiji_Affine{j},pixelSize,0);
%                 All_Dis_um = [All_Dis_um,Dis_um];
%                 All_Dis_voxel = [All_Dis_voxel,Dis_voxel];
%                 All_Dis_NonRigid_um = [All_Dis_NonRigid_um,Dis_NonRigid_um];
%                 All_Dis_NonRigid_voxel = [All_Dis_NonRigid_voxel,Dis_NonRigid_voxel];
%                 All_Dis_fiji_Affine_um = [All_Dis_fiji_Affine_um,Dis_fiji_Affine_um];
%                 All_Dis_fiji_Affine_voxel = [All_Dis_fiji_Affine_voxel,Dis_fiji_Affine_voxel];
%             end
%         end
%         hist1 = hist(All_Dis_voxel,[0.5:1:29.5]);
%         hist2 = hist(All_Dis_NonRigid_voxel,[0.5:1:29.5]);
%         hist3 = hist(All_Dis_fiji_Affine_voxel,[0.5:1:29.5]);
% %         figure,plot(hist1)
% %         hold on
% %         plot(hist2)
% %         hold on
% %         plot(hist3)
% %         legend(['trace mismatch before registration:','mean=',num2str(mean(All_Dis_voxel)),',std=',num2str(std(All_Dis_voxel))],['trace mismatch after Non-Rigid registration:','mean=',num2str(mean(All_Dis_NonRigid_voxel)),',std=',num2str(std(All_Dis_NonRigid_voxel))],['trace mismatch after fiji-Affine registration:','mean=',num2str(mean(All_Dis_fiji_Affine_voxel)),',std=',num2str(std(All_Dis_fiji_Affine_voxel))])
%     end
%     % % j = 3;
%     % % figure,PlotAM(result{1,1}.Trace.AM{j}, result{1,1}.Trace.r2{j},'b')
%     % % hold on
%     % % PlotAM(result{1,1}.Trace.AM{j}, result{1,1}.Trace.r1{j},'r')
%     % % hold on
%     % % PlotAM(result{1,1}.Trace.AM{j}, result{1,1}.Trace.r1_NR{j},'g')
%     % % legend('Target(b)','Source(r)','Non-Rigid(g)')
%     eucl_bouton = [];
%     eucl_bouton_NR = [];
%     eucl_bouton_fiji_Affine = [];
% %     figure
%     for i = 1:size(result,2)
%         d_bouton = result{1,i}.Bouton.r1-result{1,i}.Bouton.r2;
% %         plot(d_bouton(:,1),d_bouton(:,2),'.','color','r')
% %         hold on
%         d_bouton_NR = result{1,i}.Bouton.r1_NR-result{1,i}.Bouton.r2;
% %         plot(d_bouton_NR(:,1),d_bouton_NR(:,2),'.','color','g')
% %         hold on
%         d_bouton_fiji_Affine = result{1,i}.Bouton.r1_fiji_Affine-result{1,i}.Bouton.r2_fiji_Affine;
% %         plot(d_bouton_fiji_Affine(:,1),d_bouton_fiji_Affine(:,2),'.','color','b')
% %         hold on
% % %         eucl_bouton_NR1(i)=mean(sum((result{1,i}.Bouton.r2-result{1,i}.Bouton.r1_NR).^2,2).^0.5)
%         eucl_bouton = [eucl_bouton;sqrt(sum(d_bouton.^2,2))];
%         eucl_bouton_NR = [eucl_bouton_NR;sqrt(sum(d_bouton_NR.^2,2))];
%         eucl_bouton_fiji_Affine = [eucl_bouton_fiji_Affine;sqrt(sum(d_bouton_fiji_Affine.^2,2))];
%     end
%     legend('before registration','after Non-Rigid registration','after fiji-Affine registration')
%     if Calculate_Trace == 1
%         save([filename,'\Evaluation.mat'],'All_Dis_um','All_Dis_voxel','All_Dis_NonRigid_um','All_Dis_NonRigid_voxel','All_Dis_fiji_Affine_um','All_Dis_fiji_Affine_voxel','eucl_bouton','eucl_bouton_NR','eucl_bouton_fiji_Affine')
%     else
%         save([filename,'\Evaluation.mat'],'eucl_bouton','eucl_bouton_NR','eucl_bouton_fiji_Affine')
%         
%     end
%     bhist1 = hist(eucl_bouton,[0.5:1:39.5]);
%     bhist2 = hist(eucl_bouton_NR,[0.5:1:39.5]);
%     bhist3 = hist(eucl_bouton_fiji_Affine,[0.5:1:39.5]);
% %     figure
% %     plot(bhist1)
% %     hold on
% %     plot(bhist2)
% %     hold on
% %     plot(bhist3)
% %     legend(['bouton mismatch before registration:','mean=',num2str(mean(eucl_bouton)),',std=',num2str(std(eucl_bouton))],['bouton mismatch after registration:','mean=',num2str(mean(eucl_bouton_NR)),',std=',num2str(std(eucl_bouton_NR))],['bouton mismatch after registration:','mean=',num2str(mean(eucl_bouton_fiji_Affine)),',std=',num2str(std(eucl_bouton_fiji_Affine))])
%     toc
%     NR_mean(nummu) = mean(eucl_bouton_NR);
%     NR_std(nummu) = std(eucl_bouton_NR);
end
figure,errorbar([-11,-10:10],NR_mean_Trace,NR_std_Trace)
% save(['E:\Shih-Luen\Lab\Projects\RegistrationEvaluation\MatchedPoints_Non-Rigid_Seyed\','stat.mat'],'NR_mean','NR_std')
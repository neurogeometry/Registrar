filename = '..\..\RegistrationEvaluation\result1and2';
load([filename,'\Evaluation.mat'])
tic
hist1 = hist(All_Dis_voxel,[0.5:1:29.5]);
hist2 = hist(All_Dis_NonRigid_voxel,[0.5:1:29.5]);

figure
plot(hist1)
hold on
plot(hist2)
legend(['trace mismatch before registration:','mean=',num2str(mean(All_Dis_um)),',std=',num2str(std(All_Dis_um))],['trace mismatch after registration:','mean=',num2str(mean(All_Dis_NonRigid_um)),',std=',num2str(std(All_Dis_NonRigid_um))])

bhist1 = hist(eucl_bouton,[0.5:1:39.5]);
bhist2 = hist(eucl_bouton_NR,[0.5:1:39.5]);

figure
plot(bhist1)
hold on
plot(bhist2)
legend(['bouton mismatch before registration:','mean=',num2str(mean(eucl_bouton)),',std=',num2str(std(eucl_bouton))],['bouton mismatch after registration:','mean=',num2str(mean(eucl_bouton_NR)),',std=',num2str(std(eucl_bouton_NR))])
toc
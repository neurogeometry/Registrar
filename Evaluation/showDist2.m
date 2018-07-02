clear
load('E:\Shih-Luen\Lab\Projects\RegistrationEvaluation\result_NR_fiji\Evaluation.mat');
load('E:\Shih-Luen\Lab\Projects\RegistrationEvaluation\result_NR_fiji\result_NR_fiji.mat');
z = 1;
i = 1;
while z <= size(All_Dis_fiji_Affine_voxel,2)
    j = size(result{1,i}.Trace.r2,2);
    Dis_voxel{i} = All_Dis_voxel(z:z+j-1);
    Dis_NonRigid_voxel{i} = All_Dis_NonRigid_voxel(z:z+j-1);
    Dis_fiji_Affine_voxel{i} = All_Dis_fiji_Affine_voxel(z:z+j-1);
    figure
    plot(Dis_voxel{i})
    hold on
    plot(Dis_NonRigid_voxel{i})
    hold on
    plot(Dis_fiji_Affine_voxel{i})
    legend('before','after non-rigid','after fiji-Affine')
    z = z + j;
    i = i + 1;
end
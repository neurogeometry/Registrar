
mu_Names = {'0','0.00097656','0.0078125','0.0039063','0.0039063','0.015625','0.03125','0.0625','0.125','0.25','0.5','1','2','4','8','16','32','64','128','256','512','1024'};

Read_Trace = 0;
bouton_satistics = zeros(length(mu_Names),1);
for nummu = 1:size(mu_Names,2)
    filename = ['..\..\RegistrationEvaluation\MatchedPoints_Non-Rigid_Seyed\result_NR_fiji_',mu_Names{nummu}];
if Read_Trace == 1
    load([filename,'\Evaluation.mat'])
else
    load([filename,'\Evaluation.mat'])
end
mean(eucl_bouton)
std(eucl_bouton)
num2str(mean(eucl_bouton_NR))
num2str(std(eucl_bouton_NR))
mean(eucl_bouton_fiji_Affine)
std(eucl_bouton_fiji_Affine)
end
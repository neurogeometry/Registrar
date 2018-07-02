load('E:\Shih-Luen\Lab\Projects\RegistrationEvaluation\result_T0\IMs1.mat')
load('E:\Shih-Luen\Lab\Projects\RegistrationEvaluation\result_T0\result_T.mat')
figure,imshowpair(IM_source_max_P,IM_Target_NR_max_P,'Scaling','independent')
hold on
PlotAM(result{1,1}.Trace.r)
% This function performs automated tracing of neurites in a thresholded image, Im

function [AMlbl_merged,r_merged,R_merged]=AutomatedTracing(Orig,Im,Parameters,TrainingHistory_file_path)

% 1. Get Initial Trace
[AMlbl,r,R]=InitialTrace_VC(Orig,Im,Parameters);

if length(unique(AMlbl(AMlbl>0)))>1
    % 2. Generate Cluster Structure of cost components
    [AMlbl,r,R,ClustersStr]=GenerateClustersStrW(Orig,AMlbl,r,R,Parameters);
    
    % Automated Branch Merger
    temp=load(TrainingHistory_file_path);
    TrainingHistory=temp.TrainingHistory;
    [AMlbl_merged, r_merged, R_merged]=AutomatedMerger(AMlbl,r,R,ClustersStr,TrainingHistory,Parameters);
    
    % Optimize the trace
    [AMlbl_merged, r_merged, R_merged, ~]=Optimize_Trace_R(Orig,AMlbl_merged,r_merged,R_merged,0,1,0,Parameters.Optimization.PointsPerVoxel,Parameters.Optimization.MaxStepNumber,Parameters.Optimization.TraceStiffness,Parameters.Optimization.RadiusStiffness,Parameters.Optimization.StepSize,1,1);
else
    AMlbl_merged=AMlbl;
    r_merged=r;
    R_merged=R;
    disp('There is 1 tree in the stack.')
end

    

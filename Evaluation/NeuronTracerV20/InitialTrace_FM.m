% This function uses Fast Marching algorithm to calculate the initial trace required for automated tracing. 
% All branches are disconnected at branch points

function [AMlbl,r,R]=InitialTrace_FM(Orig,Im,SVr,N_steps,Max_Known_Dist,Contrast,R_min,R_step,R_max,unisotropy,Parameters,varargin)

AMlbl=[]; 
r=[]; 
R=[];

if length(varargin)==1
    myClass=varargin{1};
else
    myClass.getFlag=NaN;
end

sizeIm=size(Orig);
if length(sizeIm)==2
    sizeIm(3)=1;
end

% 1. Fast Marching
disp('Fast Marching started.')
[AMlbl,r]=FastMarching(Im,SVr,N_steps,Max_Known_Dist,Contrast,R_min,R_max,unisotropy,myClass);
disp('Fast Marching is complete.')

if myClass.getFlag()==1
    return;
end

% 2. Label branches in AM
AMlbl = LabelBranchesAM(AMlbl);
disp('Adjacency matrix is labeled.')
disp(['There are ',num2str(length(unique(AMlbl(AMlbl>0)))), ' branches in the image.'])
%AMlbl = LabelTreesAM(AMlbl);???

if myClass.getFlag()==1
    return;
end

% 3. Fixing the tips of terminal branches
[AMlbl,r]=ExtendTips(Orig,AMlbl,r,ceil(Max_Known_Dist*0.75),unisotropy);
[AMlbl,r]=TrimTips(AMlbl,r,Orig,Max_Known_Dist,R_min,max(0.2,(R_max-R_min)/10),R_max);
disp('Tips of terminal branches are fixed.')

% % 4. Eliminate Short Vertical Intermediate Branches
% for i=2:2:Parameters.FastMarching.MinIntermediateBranchLength
%     if myClass.getFlag()==1
%         return;
%     end
%     [AMlbl,r] = Eliminate_Vertical_Intermediate_Branches(AMlbl,r,i);
% end
% disp([num2str(length(unique(AMlbl(AMlbl>0)))), ' branches remain after short vertical intermediate branch elimination.'])

% % 5. Eliminate Short Terminal Branches
% for thr=1:1:Parameters.FastMarching.MinTerminalBranchLength
%     if myClass.getFlag()==1
%         return;
%     end
%     [AMlbl r ~] = Eliminate_Terminal_Branches(AMlbl,r,thr,1,0);     
% end
% disp([num2str(length(unique(AMlbl(AMlbl>0)))), ' branches remain after short terminal branch elimination.'])

if myClass.getFlag()==1
    return;
end
    
% 6. Eliminate Small Loops
disp('Eliminating small loops.')
[AMlbl, r, ~, ~]=Reduce_Small_Loops(AMlbl,r,Parameters.FastMarching.MinLoopSize);
disp([num2str(length(unique(AMlbl(AMlbl>0)))), ' branches remain after small loop elimination.'])
   
% 7. Eliminate Short Terminal Branches
for thr=2:2:Parameters.FastMarching.MinTerminalBranchLength
    if myClass.getFlag()==1
        return;
    end
    [AMlbl r ~] = Eliminate_Terminal_Branches(AMlbl,r,thr,1,0);    
end
disp([num2str(length(unique(AMlbl(AMlbl>0)))), ' branches remain after short terminal branch elimination.'])

% 8. Reduce or Eliminate Short Intermediate Branches
for i=2:2:Parameters.FastMarching.MinIntermediateBranchLength
    if myClass.getFlag()==1
        return;
    end
    [AMlbl,r] = Reduce_Short_Intermediate_Branches(AMlbl,r,i);
    %[AMlbl,r] = Eliminate_Short_Intermediate_Branches(AMlbl,r,i);    
end
disp([num2str(length(unique(AMlbl(AMlbl>0)))), ' branches remain after short intermediate branch reduction.'])

% 9. Cut the trace near the xz and yz faces of the stack
edge_points=(r(:,1)<=Parameters.FastMarching.TrimTrace | r(:,1)>sizeIm(1)-Parameters.FastMarching.TrimTrace | r(:,2)<=Parameters.FastMarching.TrimTrace | r(:,2)>sizeIm(2)-Parameters.FastMarching.TrimTrace);
if nnz(edge_points)>0
    AMlbl(edge_points,:)=[];
    AMlbl(:,edge_points)=[];
    r(edge_points,:)=[];
    AMlbl = LabelBranchesAM(AMlbl>0);
end

if myClass.getFlag()==1
    return;
end
    
% 10. Optimize the trace
[R ~] = Find_R_fast(Orig,r,R_min,max(0.2,(R_max-R_min)/10),R_max);
[AMlbl, r, R, ~]=Optimize_Trace_R(Orig,AMlbl,r,R,0,0,0,Parameters.Optimization.PointsPerVoxel,Parameters.Optimization.MaxStepNumber,Parameters.Optimization.TraceStiffness,Parameters.Optimization.RadiusStiffness,Parameters.Optimization.StepSize,1,1,myClass);
% AMlbl is labeled for trees after this step
 
disp('Initial trace is created.')


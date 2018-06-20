% This function uses Voxel coding algorithm to calculate the initial trace required for automated tracing. 
% All branches are disconnected at branch points

function [AMlbl,r,R]=InitialTrace_VC(Orig,Im,Parameters,varargin)

AMlbl=[]; 
r=[]; 
R=[];

if length(varargin)==1
    myClass=varargin{1};
else
    myClass.getFlag=NaN;
end

disp(['There are ',num2str(nnz(Im)), ' voxels in the thresholded image.'])

sizeIm=size(Im);
if length(sizeIm)==2
    sizeIm(3)=1;
end

% 1. Voxel Coding
disp('Voxel Coding started.')
Start_X=[]; Start_Y=[]; Start_Z=[];
W = VoxelCoding(Im,Start_X,Start_Y,Start_Z,1,myClass);
disp('Voxel Coding is complete.')

if myClass.getFlag()==1
    return;
end

% 2. Create sparse AM from W
[AM,r] = W2AM(Im,W,[]);
[AM r] = Remove_Zerolength_Segments(AM,r);
clear W

disp('Adjacency matrix is created.')

if myClass.getFlag()==1
    return;
end

% 3. Label branches in AM
AMlbl = LabelBranchesAM(AM);
disp('Adjacency matrix is labeled.')
disp(['There are ',num2str(length(unique(AMlbl(AMlbl>0)))), ' branches in the image.'])

% 4. Eliminate Short Vertical Intermediate Branches
for i=2:2:Parameters.VoxelCoding.MinIntermediateBranchLength
    if myClass.getFlag()==1
        return;
    end
    [AMlbl,r] = Eliminate_Vertical_Intermediate_Branches(AMlbl,r,i);
end
disp([num2str(length(unique(AMlbl(AMlbl>0)))), ' branches remain after short vertical intermediate branch elimination.'])

% 5. Eliminate Short Terminal Branches
for thr=1:1:5
    if myClass.getFlag()==1
        return;
    end
    [AMlbl r ~] = Eliminate_Terminal_Branches(AMlbl,r,thr,1,0);     
end
disp([num2str(length(unique(AMlbl(AMlbl>0)))), ' branches remain after short terminal branch elimination.'])

if myClass.getFlag()==1
    return;
end
    
% 6. Eliminate Small Loops
disp('Eliminating small loops.')
[AMlbl, r, ~, ~]=Reduce_Small_Loops(AMlbl,r,Parameters.VoxelCoding.MinLoopSize);
disp([num2str(length(unique(AMlbl(AMlbl>0)))), ' branches remain after small loop elimination.'])
   
% 7. Eliminate Short Terminal Branches
for thr=6:1:Parameters.VoxelCoding.MinTerminalBranchLength
    if myClass.getFlag()==1
        return;
    end
    [AMlbl r ~] = Eliminate_Terminal_Branches(AMlbl,r,thr,1,0);    
end
disp([num2str(length(unique(AMlbl(AMlbl>0)))), ' branches remain after short terminal branch elimination.'])

% 8. Reduce or Eliminate Short Intermediate Branches
for i=2:2:Parameters.VoxelCoding.MinIntermediateBranchLength
    if myClass.getFlag()==1
        return;
    end
    [AMlbl,r] = Reduce_Short_Intermediate_Branches(AMlbl,r,i);
    %[AMlbl,r] = Eliminate_Short_Intermediate_Branches(AMlbl,r,i);    
end
disp([num2str(length(unique(AMlbl(AMlbl>0)))), ' branches remain after short intermediate branch reduction.'])

% 9. Cut the trace near the xz and yz faces of the stack
edge_points=(r(:,1)<=Parameters.VoxelCoding.TrimTrace | r(:,1)>sizeIm(1)-Parameters.VoxelCoding.TrimTrace | r(:,2)<=Parameters.VoxelCoding.TrimTrace | r(:,2)>sizeIm(2)-Parameters.VoxelCoding.TrimTrace);
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
[AMlbl, r, R, ~]=Optimize_Trace_R(Orig,AMlbl,r,[],Parameters.Optimization.TypicalRadius,1,1,Parameters.Optimization.PointsPerVoxel,Parameters.Optimization.MaxStepNumber,Parameters.Optimization.TraceStiffness,Parameters.Optimization.RadiusStiffness,Parameters.Optimization.StepSize,1,1,myClass);
% AMlbl is labeled for trees after this step

disp('Initial trace is created.')


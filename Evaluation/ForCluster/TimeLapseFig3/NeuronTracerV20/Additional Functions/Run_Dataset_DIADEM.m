function Run_Dataset_DIADEM(dataset,stack,start_point)

data_format='Full'; 
relative_thr=0.1;
reduct_type='xyz'; % 'NA','xy','xyz'
reduct_factor=2;
reduct_method='Median';


% 1. Import thesholded image stack
switch dataset
    case 'L6'
        pth0='C:\Armen\DIADEM\Data\L6\Image Stacks\';
        pth=[pth0,stack,'\'];
        
        Start_X=1; Start_Y=1; Start_Z=1;
        
        [Orig,sizeOrig,classOrig]=ImportStack(pth,data_format,relative_thr,reduct_type,reduct_factor,reduct_method);
        
        S=1.5; % Center Surround Filter size 1,2
        thr_low=25; % Threshold 25,8
        Region_size_thr=200; % Eliminate small regions
        
        Parameters.W=1; % Propagate the wave into all or one region: W=1,0
        Parameters.ShortTerminalBranchLength=10; % Remove terminal branches that are shorter than: ShortTerminalBranchLength
        Parameters.LongBranchLength=20; % Eliminate loops that are smaller than: LongBranchLength
        Parameters.ShortIntermediateBranchLength=10; % Remove intermediate branches that are shorter than: ShortIntermediateBranchLength
        Parameters.trim=5; % Cut the trace near xy the faces of the stack
        Parameters.pointsperum=0.25; % Number of points per um in the trace: pointsperum
        Parameters.Nstep=20; % Number of trace optimization steps: Nstep
        Parameters.alpha=0.2; % Stiffness of the trace: alpha
        Parameters.betta=0.1; % Trace optimization step size: betta
        Parameters.sig_ips=1; % Intermediate-point sigma
        Parameters.sig_bps=3; % Branch-point sigma
        Parameters.sig_tps=3; % Terminal-point sigma
        Parameters.Disconnect=0; % Disconnect all branches at branch points: Disconnect=1,0
        Parameters.MaxDistance=17; % Consider merging branches that are closer than: MaxDistance
        Parameters.stepsback=5; % Number of steps taken to determine branch tip directions: stepsback
        Parameters.Dist_thr=17; % Minimal distance between branch tip clusters: Dist_thr
        % Merger cost function parameters:
        Parameters.Cos2=-1; Parameters.d2=0.1; Parameters.c2=1; Parameters.m2=-1; Parameters.f2=1.5;
        Parameters.Cos3=-0.5; Parameters.d3=0.2; Parameters.c3=1; Parameters.m3=-1; Parameters.f3=1;
        Parameters.Cos4=0; Parameters.d4=100; Parameters.c4=1; Parameters.m4=0; Parameters.f4=0;
        Parameters.g=1.5;
        Parameters.SmallTreeLength=0;
        
    case 'Olfactory'
        pth0='C:\Armen\DIADEM\Data\Olfactory_Projection_Fibers\Image Stacks\';
        pth=[pth0,stack,'\'];
        
        R_start=[119.76	215.98	39
                118.64	181.34	55
                64.56	364.47	4];
        
        Start_X=R_start(start_point,1);
        Start_Y=R_start(start_point,2);
        Start_Z=R_start(start_point,3);
         
        [Orig,sizeOrig,classOrig]=ImportStack(pth,data_format,relative_thr,reduct_type,reduct_factor,reduct_method);
        
        if strcmp(stack,'OP_A') || strcmp(stack,'OP_B') || strcmp(stack,'OP_C')
            S=[1,1.5,2]; % Center Surround Filter size
            thr_low=15; % Threshold
            Region_size_thr=200; % Eliminate small regions
            
            Parameters.W=1; % Propagate the wave into all or one region: W=1,0
            Parameters.ShortTerminalBranchLength=6; % Remove terminal branches that are shorter than: ShortTerminalBranchLength
            Parameters.LongBranchLength=15; % Eliminate loops that are smaller than: LongBranchLength
            Parameters.ShortIntermediateBranchLength=0; % Remove intermediate branches that are shorter than: ShortIntermediateBranchLength
            Parameters.trim=0; % Cut the trace near the xy faces of the stack
            Parameters.pointsperum=0.25; % Number of points per um in the trace: pointsperum
            Parameters.Nstep=20; % Number of trace optimization steps: Nstep
            Parameters.alpha=0.1; % Stiffness of the trace: alpha
            Parameters.betta=0.1; % Trace optimization step size: betta
            Parameters.sig_ips=1; % Intermediate-point sigma
            Parameters.sig_bps=3; % Branch-point sigma
            Parameters.sig_tps=3; % Terminal-point sigma
            Parameters.Disconnect=1; % Disconnect all branches at branch points: Disconnect=1,0
            Parameters.MaxDistance=8; % Consider merging branches that are closer than: MaxDistance
            Parameters.stepsback=3; % Number of steps taken to determine branch tip directions: stepsback
            Parameters.Dist_thr=8; % Minimal distance between branch tip clusters: Dist_thr
            % Merger cost function parameters:
            Parameters.Cos2=-1; Parameters.d2=0.1; Parameters.c2=2; Parameters.m2=-1; Parameters.f2=10;
            Parameters.Cos3=-0.5; Parameters.d3=0.1; Parameters.c3=0.5; Parameters.m3=-1; Parameters.f3=0;
            Parameters.Cos4=0; Parameters.d4=0.1; Parameters.c4=0; Parameters.m4=-3; Parameters.f4=0;
            Parameters.g=1.5;
            Parameters.SmallTreeLength=0;
        elseif strcmp(stack,'OP_4')
            S=1; % Center Surround Filter size
            thr_low=20; % Threshold
            Region_size_thr=200; % Eliminate small regions
            
            Parameters.W=1; % Propagate the wave into all or one region: W=1,0
            Parameters.ShortTerminalBranchLength=8; % Remove terminal branches that are shorter than: ShortTerminalBranchLength
            Parameters.LongBranchLength=15; % Eliminate loops that are smaller than: LongBranchLength
            Parameters.ShortIntermediateBranchLength=0; % Remove intermediate branches that are shorter than: ShortIntermediateBranchLength
            Parameters.trim=2; % Cut the trace near the xy faces of the stack
            Parameters.pointsperum=0.25; % Number of points per um in the trace: pointsperum
            Parameters.Nstep=20; % Number of trace optimization steps: Nstep
            Parameters.alpha=0.1; % Stiffness of the trace: alpha
            Parameters.betta=0.1; % Trace optimization step size: betta
            Parameters.sig_ips=1; % Intermediate-point sigma
            Parameters.sig_bps=3; % Branch-point sigma
            Parameters.sig_tps=3; % Terminal-point sigma
            Parameters.Disconnect=1; % Disconnect all branches at branch points: Disconnect=1,0
            Parameters.MaxDistance=5; % Consider merging branches that are closer than: MaxDistance
            Parameters.stepsback=20; % Number of steps taken to determine branch tip directions: stepsback
            Parameters.Dist_thr=3; % Minimal distance between branch tip clusters: Dist_thr
            % Merger cost function parameters:
            Parameters.Cos2=-1; Parameters.d2=0.1; Parameters.c2=2; Parameters.m2=-1; Parameters.f2=0;
            Parameters.Cos3=-0.5; Parameters.d3=0.1; Parameters.c3=0.5; Parameters.m3=-1; Parameters.f3=0;
            Parameters.Cos4=0; Parameters.d4=0.1; Parameters.c4=1; Parameters.m4=0; Parameters.f4=0;
            Parameters.g=2;
            Parameters.SmallTreeLength=0;
        elseif strcmp(stack,'OP_5')
            S=0; % Center Surround Filter size
            thr_low=16; % Threshold
            Region_size_thr=200; % Eliminate small regions
            
            Parameters.W=1; % Propagate the wave into all or one region: W=1,0
            Parameters.ShortTerminalBranchLength=10; % Remove terminal branches that are shorter than: ShortTerminalBranchLength
            Parameters.LongBranchLength=15; % Eliminate loops that are smaller than: LongBranchLength
            Parameters.ShortIntermediateBranchLength=0; % Remove intermediate branches that are shorter than: ShortIntermediateBranchLength
            Parameters.trim=2; % Cut the trace near the xy faces of the stack
            Parameters.pointsperum=0.25; % Number of points per um in the trace: pointsperum
            Parameters.Nstep=20; % Number of trace optimization steps: Nstep
            Parameters.alpha=0.1; % Stiffness of the trace: alpha
            Parameters.betta=0.1; % Trace optimization step size: betta
            Parameters.sig_ips=1; % Intermediate-point sigma
            Parameters.sig_bps=3; % Branch-point sigma
            Parameters.sig_tps=3; % Terminal-point sigma
            Parameters.Disconnect=0; % Disconnect all branches at branch points: Disconnect=1,0
            Parameters.MaxDistance=5; % Consider merging branches that are closer than: MaxDistance
            Parameters.stepsback=5; % Number of steps taken to determine branch tip directions: stepsback
            Parameters.Dist_thr=5; % Minimal distance between branch tip clusters: Dist_thr
            % Merger cost function parameters:
            Parameters.Cos2=-1; Parameters.d2=0.1; Parameters.c2=2; Parameters.m2=-1; Parameters.f2=0;
            Parameters.Cos3=-0.5; Parameters.d3=0.1; Parameters.c3=0.5; Parameters.m3=-1; Parameters.f3=0;
            Parameters.Cos4=0; Parameters.d4=0.1; Parameters.c4=1; Parameters.m4=0; Parameters.f4=0;
            Parameters.g=2;
            Parameters.SmallTreeLength=0;
        elseif strcmp(stack,'OP_6') % !!! perfect parameters
            S=0; % Center Surround Filter size
            thr_low=20; % Threshold
            Region_size_thr=200; % Eliminate small regions
            
            Parameters.W=1; % Propagate the wave into all or one region: W=1,0
            Parameters.ShortTerminalBranchLength=6; % Remove terminal branches that are shorter than: ShortTerminalBranchLength
            Parameters.LongBranchLength=15; % Eliminate loops that are smaller than: LongBranchLength
            Parameters.ShortIntermediateBranchLength=0; % Remove intermediate branches that are shorter than: ShortIntermediateBranchLength
            Parameters.trim=2; % Cut the trace near the xy faces of the stack
            Parameters.pointsperum=0.25; % Number of points per um in the trace: pointsperum
            Parameters.Nstep=20; % Number of trace optimization steps: Nstep
            Parameters.alpha=0.1; % Stiffness of the trace: alpha
            Parameters.betta=0.1; % Trace optimization step size: betta
            Parameters.sig_ips=1; % Intermediate-point sigma
            Parameters.sig_bps=3; % Branch-point sigma
            Parameters.sig_tps=3; % Terminal-point sigma
            Parameters.Disconnect=0; % Disconnect all branches at branch points: Disconnect=1,0
            Parameters.MaxDistance=5; % Consider merging branches that are closer than: MaxDistance
            Parameters.stepsback=5; % Number of steps taken to determine branch tip directions: stepsback
            Parameters.Dist_thr=5; % Minimal distance between branch tip clusters: Dist_thr
            % Merger cost function parameters:
            Parameters.Cos2=-1; Parameters.d2=0.1; Parameters.c2=2; Parameters.m2=-1; Parameters.f2=0;
            Parameters.Cos3=-0.5; Parameters.d3=0.1; Parameters.c3=0.5; Parameters.m3=-1; Parameters.f3=0;
            Parameters.Cos4=0; Parameters.d4=0.1; Parameters.c4=1; Parameters.m4=0; Parameters.f4=0;
            Parameters.g=2;
            Parameters.SmallTreeLength=0;
        end
        
    case 'Cerebellar'
        pth0='C:\Armen\DIADEM\Data\Cerebellar Climbing Fibers\Image Stacks\'; 
        pth=[pth0,stack,'\'];
        
        Start_X=1138;  
        Start_Y=773;
        Start_Z=11;
        
       [Orig,sizeOrig,classOrig]=ImportStack(pth,data_format,relative_thr,reduct_type,reduct_factor,reduct_method);
        
        Start_X=(Start_X/reduct_factor);
        Start_Y=(Start_Y/reduct_factor);
        Start_Z=(Start_Z/reduct_factor);      
       
        S=[2,3,4]; % Center Surround Filter size
        thr_low=20; % Threshold
        Region_size_thr=200; % Eliminate small regions
        
        Parameters.W=1; % Propagate the wave into all or one region: W=1,0
        Parameters.ShortTerminalBranchLength=8; % Remove terminal branches that are shorter than: ShortTerminalBranchLength
        Parameters.LongBranchLength=20; % Eliminate loops that are smaller than: LongBranchLength
        Parameters.ShortIntermediateBranchLength=0; % Remove intermediate branches that are shorter than: ShortIntermediateBranchLength
        Parameters.trim=2; % Cut the trace near the xy faces of the stack
        Parameters.pointsperum=0.25; % Number of points per um in the trace: pointsperum
        Parameters.Nstep=20; % Number of trace optimization steps: Nstep
        Parameters.alpha=0.1; % Stiffness of the trace: alpha
        Parameters.betta=0.1; % Trace optimization step size: betta
        Parameters.sig_ips=1; % Intermediate-point sigma
        Parameters.sig_bps=3; % Branch-point sigma
        Parameters.sig_tps=3; % Terminal-point sigma
        Parameters.Disconnect=0; % Disconnect all branches at branch points: Disconnect=1,0
        Parameters.MaxDistance=7; % Consider merging branches that are closer than: MaxDistance
        Parameters.stepsback=5; % Number of steps taken to determine branch tip directions: stepsback
        Parameters.Dist_thr=7; % Minimal distance between branch tip clusters: Dist_thr
        % Merger cost function parameters:
        Parameters.Cos2=-1; Parameters.d2=0.1; Parameters.c2=2; Parameters.m2=-1; Parameters.f2=0;
        Parameters.Cos3=-0.5; Parameters.d3=0.1; Parameters.c3=0.5; Parameters.m3=-1; Parameters.f3=0;
        Parameters.Cos4=0; Parameters.d4=100; Parameters.c4=1; Parameters.m4=0; Parameters.f4=0;
        Parameters.g=5;
        Parameters.SmallTreeLength=0;
        
    case 'CA3'
        pth0='C:\Armen\DIADEM\Data\CA3\Image Stacks\';
        pth=[pth0,stack,'\'];
        
        [Orig,sizeOrig,classOrig]=ImportStack(pth,data_format,relative_thr,reduct_type,reduct_factor,reduct_method);
        
        Start_X=1; Start_Y=1; Start_Z=1;     
       
        S=2; % Center Surround Filter size
        thr_low=50; % Threshold
        Region_size_thr=200; % Eliminate small regions
        
        Parameters.W=1; % Propagate the wave into all or one region: W=1,0
        Parameters.ShortTerminalBranchLength=10; % Remove terminal branches that are shorter than: ShortTerminalBranchLength
        Parameters.LongBranchLength=15; % Eliminate loops that are smaller than: LongBranchLength
        Parameters.ShortIntermediateBranchLength=1; % Remove intermediate branches that are shorter than: ShortIntermediateBranchLength
        Parameters.trim=2; % Cut the trace near the xy faces of the stack
        Parameters.pointsperum=0.25; % Number of points per um in the trace: pointsperum
        Parameters.Nstep=20; % Number of trace optimization steps: Nstep
        Parameters.alpha=0.2; % Stiffness of the trace: alpha
        Parameters.betta=0.1; % Trace optimization step size: betta
        Parameters.sig_ips=1; % Intermediate-point sigma
        Parameters.sig_bps=3; % Branch-point sigma
        Parameters.sig_tps=3; % Terminal-point sigma
        Parameters.Disconnect=0; % Disconnect all branches at branch points: Disconnect=1,0
        Parameters.MaxDistance=5; % Consider merging branches that are closer than: MaxDistance
        Parameters.stepsback=50; % Number of steps taken to determine branch tip directions: stepsback
        Parameters.Dist_thr=5; % Minimal distance between branch tip clusters: Dist_thr
        % Merger cost function parameters:
        Parameters.Cos2=-1; Parameters.d2=0.1; Parameters.c2=2; Parameters.m2=-1; Parameters.f2=0;
        Parameters.Cos3=-0.5; Parameters.d3=0.1; Parameters.c3=0.5; Parameters.m3=-1; Parameters.f3=0;
        Parameters.Cos4=0; Parameters.d4=100; Parameters.c4=1; Parameters.m4=0; Parameters.f4=0;
        Parameters.g=2;
        Parameters.SmallTreeLength=0;

    case 'Neuromuscular'
        pth0='C:\Armen\DIADEM\Data\Neuromuscular\Image Stacks\';
        pth=[pth0,stack,'\'];
        
        [Orig,sizeOrig,classOrig]=ImportStack(pth,data_format,relative_thr,reduct_type,reduct_factor,reduct_method);
        Orig=255-Orig;
        
        Start_X=1; Start_Y=1; Start_Z=1;
        
        S=1; % Center Surround Filter size
        thr_low=10; % Threshold
        Region_size_thr=100; % Eliminate small regions
        
        Parameters.W=1; % Propagate the wave into all or one region: W=1,0
        Parameters.ShortTerminalBranchLength=5; % Remove terminal branches that are shorter than: ShortTerminalBranchLength
        Parameters.LongBranchLength=20; % Eliminate loops that are smaller than: LongBranchLength
        Parameters.ShortIntermediateBranchLength=7; % Remove intermediate branches that are shorter than: ShortIntermediateBranchLength
        Parameters.trim=5; % Cut the trace near the xy faces of the stack
        Parameters.pointsperum=0.25; % Number of points per um in the trace: pointsperum
        Parameters.Nstep=20; % Number of trace optimization steps: Nstep
        Parameters.alpha=0.5; % Stiffness of the trace: alpha
        Parameters.betta=0.1; % Trace optimization step size: betta
        Parameters.sig_ips=1; % Intermediate-point sigma
        Parameters.sig_bps=3; % Branch-point sigma
        Parameters.sig_tps=3; % Terminal-point sigma
        Parameters.Disconnect=0; % Disconnect all branches at branch points: Disconnect=1,0
        Parameters.MaxDistance=12; % Consider merging branches that are closer than: MaxDistance
        Parameters.stepsback=5; % Number of steps taken to determine branch tip directions: stepsback
        Parameters.Dist_thr=15; % Minimal distance between branch tip clusters: Dist_thr
        % Merger cost function parameters:
        Parameters.Cos2=-1; Parameters.d2=0.1; Parameters.c2=2; Parameters.m2=-0; Parameters.f2=10;
        Parameters.Cos3=-0.5; Parameters.d3=0.2; Parameters.c3=2; Parameters.m3=-0; Parameters.f3=10;
        Parameters.Cos4=0; Parameters.d4=100; Parameters.c4=1; Parameters.m4=0; Parameters.f4=0;
        Parameters.g=2;
        Parameters.SmallTreeLength=30;
   
        
    case 'V1'
        pth0='C:\Armen\DIADEM\Data\Visual Cortex\Image Stacks\Section 1\';
        %pth0='C:\Armen\DIADEM\Data\Visual Cortex\Image Stacks\Section 2\';
        %pth0='C:\Armen\DIADEM\Data\Visual Cortex\Image Stacks\Section 3\';
        %pth0='C:\Armen\DIADEM\Data\Visual Cortex\Image Stacks\Section 4\';
        pth=[pth0,stack,'\'];
        
        [Orig,sizeOrig,classOrig]=ImportStack(pth,data_format,relative_thr,reduct_type,reduct_factor,reduct_method);
                
        Start_X=1; Start_Y=1; Start_Z=1;
        
        S=1; % Center Surround Filter size
        thr_low=30; % Threshold
        Region_size_thr=200; % Eliminate small regions
        
        Parameters.W=1; % Propagate the wave into all or one region: W=1,0
        Parameters.ShortTerminalBranchLength=15; % Remove terminal branches that are shorter than: ShortTerminalBranchLength
        Parameters.LongBranchLength=20; % Eliminate loops that are smaller than: LongBranchLength
        Parameters.ShortIntermediateBranchLength=5; % Remove intermediate branches that are shorter than: ShortIntermediateBranchLength
        Parameters.trim=5; % Cut the trace near the xy faces of the stack
        Parameters.pointsperum=0.25; % Number of points per um in the trace: pointsperum
        Parameters.Nstep=20; % Number of trace optimization steps: Nstep
        Parameters.alpha=0.5; % Stiffness of the trace: alpha
        Parameters.betta=0.1; % Trace optimization step size: betta
        Parameters.sig_ips=1; % Intermediate-point sigma
        Parameters.sig_bps=3; % Branch-point sigma
        Parameters.sig_tps=3; % Terminal-point sigma
        Parameters.Disconnect=0; % Disconnect all branches at branch points: Disconnect=1,0
        Parameters.MaxDistance=12; % Consider merging branches that are closer than: MaxDistance
        Parameters.stepsback=10; % Number of steps taken to determine branch tip directions: stepsback
        Parameters.Dist_thr=15; % Minimal distance between branch tip clusters: Dist_thr
        % Merger cost function parameters:
        Parameters.Cos2=-1; Parameters.d2=0.1; Parameters.c2=2; Parameters.m2=-0; Parameters.f2=10;
        Parameters.Cos3=-0.5; Parameters.d3=0.2; Parameters.c3=2; Parameters.m3=-0; Parameters.f3=10;
        Parameters.Cos4=0; Parameters.d4=100; Parameters.c4=1; Parameters.m4=0; Parameters.f4=0;
        Parameters.g=2;
        Parameters.SmallTreeLength=0;
end

if strcmp(data_format,'Full') 
    %Im = ReduceBackground(Orig,10);
    Im=CSF_cos(Orig,S); 
    Im=Threshold(Im,thr_low); 
    Im=EliminateSmallReg(Im,Region_size_thr);
    [AMlbl_merged,r_merged]=AutomatedTracing(Im,Start_X,Start_Y,Start_Z,Parameters);
end


% 3. Ploting the trace on a maximum projection image
figure
imshow(max(Orig,[],3),[0 max(Orig(:))]) 
hold on 
PlotAM(AMlbl_merged, r_merged)


% 4. Exporting the trace
save_pth=[pth0,'Results\'];
[swc_trace_1,swc_trace_all] = AM2swc(AMlbl_merged,r_merged,reduct_type,reduct_factor); 
dlmwrite([save_pth,stack,'_1.swc'], swc_trace_1, 'delimiter', ' ')
dlmwrite([save_pth,stack,'_all.swc'], swc_trace_all, 'delimiter', ' ')

IM=Orig;
AM=AMlbl_merged;
R=r_merged;
Red_Type=reduct_type;
Red_Fac=reduct_factor;

% 5. Save AMlbl_merged,r_merged,reduct_type,reduct_factor
save([save_pth,stack,'_AM-R.mat'],'IM','AM','R','Red_Type','Red_Fac')


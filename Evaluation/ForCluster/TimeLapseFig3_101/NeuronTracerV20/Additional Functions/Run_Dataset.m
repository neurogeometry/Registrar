function Run_Dataset(dataset,stack,start_point)

data_format='Full'; 
relative_thr=0.1;
reduct_type='NA'; % 'NA','xy','xyz'
reduct_factor=1;
reduct_method='Max';


% 1. Import thesholded image stack
switch dataset
    case 'L6'
        pth='C:\Armen\DIADEM\Data\L6\Image Stacks\';
        if ~strcmp(stack,'mosaic')
            %load([pth,'\Positions.mat'])
%             Start_X=1/reduct_factor; 
%             Start_Y=406/reduct_factor; 
%             Start_Z=44/reduct_factor;
            Start_X=1; Start_Y=1; Start_Z=1;

            pth=[pth,stack,'\'];
            [Orig,sizeOrig,classOrig]=ImportStack(pth,data_format,relative_thr,reduct_type,reduct_factor,reduct_method);
        elseif strcmp(stack,'mosaic')
            load([pth,'\Mosaic.mat']); % load([pth,'\Mosaic_Sparse_40.mat'])
            %Start_X=SP(start_point,1); Start_Y=SP(start_point,2); Start_Z=SP(start_point,3);
            
            %Orig=ReduceMatImage(pth,n,'Max');
            load([pth,'\Positions.mat'])
            Start_X=(SP(start_point,1)/reduct_factor);
            Start_Y=(SP(start_point,2)/reduct_factor);
            Start_Z=(SP(start_point,3)/reduct_factor);
        end
        
        S=1; % Center Surround Filter size 1,2
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
%         Parameters.Cos2=-1; Parameters.d2=0; Parameters.c2=1; Parameters.m2=-1; Parameters.f2=0;
%         Parameters.Cos3=-0.5; Parameters.d3=0.1; Parameters.c3=1; Parameters.m3=-1; Parameters.f3=0;
%         Parameters.Cos4=0; Parameters.d4=100; Parameters.c4=1; Parameters.m4=0; Parameters.f4=0;
%         Parameters.g=2;
        
%         Parameters.Cos2=-1; Parameters.d2=0.2331; Parameters.c2=1.2448; Parameters.m2=-1.1499; Parameters.f2=0.2180;
%         Parameters.Cos3=-0.5; Parameters.d3=0.3193; Parameters.c3=1.2397; Parameters.m3=0; Parameters.f3=0.9809;
%         Parameters.Cos4=0; Parameters.d4=100; Parameters.c4=1; Parameters.m4=0; Parameters.f4=0;

        Parameters.Cos2=-1; Parameters.d2=0.1; Parameters.c2=1; Parameters.m2=-1; Parameters.f2=1.5;
        Parameters.Cos3=-0.5; Parameters.d3=0.2; Parameters.c3=1; Parameters.m3=-1; Parameters.f3=1;
        Parameters.Cos4=0; Parameters.d4=100; Parameters.c4=1; Parameters.m4=0; Parameters.f4=0;
        Parameters.g=2;
        Parameters.SmallTreeLength=0;
        
    case 'Olfactory'
        %pth='C:\Armen\DIADEM\Data\Olfactory_Projection_Fibers\Image Stacks\Training Image Stack\';
        pth='C:\Armen\DIADEM\Data\Olfactory_Projection_Fibers\Image Stacks\Qualifier Image Stack\';
        Start_X=[30.979,0.72501,93.742,128.2,185.7,15.074];
        Start_Y=[429.04,391.08,179,504.37,264.02,412.01];
        Start_Z=[0,25,38,0.303,33,10];
        
        pth=[pth,stack,'\'];
        [Orig,sizeOrig,classOrig]=ImportStack(pth,data_format,relative_thr,reduct_type,reduct_factor,reduct_method);
        Start_X=Start_X(start_point);
        Start_Y=Start_Y(start_point);
        Start_Z=Start_Z(start_point);
        
        if strcmp(stack,'OP_1') || strcmp(stack,'OP_2')
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
        elseif strcmp(stack,'OP_3')
%             S=0; % Center Surround Filter size
%             thr_low=20; % Threshold
%             Region_size_thr=200; % Eliminate small regions
%             
%             Parameters.W=1; % Propagate the wave into all or one region: W=1,0
%             Parameters.ShortTerminalBranchLength=10; % Remove terminal branches that are shorter than: ShortTerminalBranchLength
%             Parameters.LongBranchLength=15; % Eliminate loops that are smaller than: LongBranchLength
%             Parameters.ShortIntermediateBranchLength=0; % Remove intermediate branches that are shorter than: ShortIntermediateBranchLength
%             Parameters.trim=2; % Cut the trace near the xy faces of the stack
%             Parameters.pointsperum=1; % Number of points per um in the trace: pointsperum
%             Parameters.Nstep=20; % Number of trace optimization steps: Nstep
%             Parameters.alpha=1; % Stiffness of the trace: alpha
%             Parameters.betta=0.05; % Trace optimization step size: betta
%             Parameters.sig_ips=1; % Intermediate-point sigma
%             Parameters.sig_bps=3; % Branch-point sigma
%             Parameters.sig_tps=3; % Terminal-point sigma
%             Parameters.Disconnect=0; % Disconnect all branches at branch points: Disconnect=1,0
%             Parameters.MaxDistance=5; % Consider merging branches that are closer than: MaxDistance
%             Parameters.stepsback=20; % Number of steps taken to determine branch tip directions: stepsback
%             Parameters.Dist_thr=5; % Minimal distance between branch tip clusters: Dist_thr
%             % Merger cost function parameters:
%             Parameters.Cos2=-1; Parameters.d2=0.1; Parameters.c2=2; Parameters.m2=-1; Parameters.f2=0;
%             Parameters.Cos3=-0.5; Parameters.d3=0.1; Parameters.c3=0.5; Parameters.m3=-1; Parameters.f3=0;
%             Parameters.Cos4=0; Parameters.d4=0.1; Parameters.c4=1; Parameters.m4=0; Parameters.f4=0;
%             Parameters.g=2;
%             Parameters.SmallTreeLength=0;
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
        %pth='C:\Armen\DIADEM\Data\Cerebellar Climbing Fibers\Image Stacks\Training Image Stack\';
        pth='C:\Armen\DIADEM\Data\Cerebellar Climbing Fibers\Image Stacks\Qualifier Image Stack\';
        Start_X=[4882,818];
        Start_Y=[1797,332];
        Start_Z=[19,6];
        
        pth=[pth,stack,'\'];
        [Orig,sizeOrig,classOrig]=ImportStack(pth,data_format,relative_thr,reduct_type,reduct_factor,reduct_method);
        
        Start_X=(Start_X(start_point)/reduct_factor);
        Start_Y=(Start_Y(start_point)/reduct_factor);
        Start_Z=(Start_Z(start_point)/reduct_factor);      
       
        S=2; % Center Surround Filter size
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
        
        %pth='C:\Armen\DIADEM\Data\CA3\Training Image Stack\';
        pth='C:\Armen\DIADEM\Data\CA3\Qualifier Image Stacks\';

        pth=[pth,stack,'\'];
        [Orig,sizeOrig,classOrig]=ImportStack(pth,data_format,relative_thr,reduct_type,reduct_factor,reduct_method);
        
        if strcmp(stack,'Section_01')
            Start_X=10; Start_Y=10; Start_Z=10;
        elseif strcmp(stack,'Section_02')
            Start_X=[4183,4503,4063,3907,3889,3613,3357,1912,2466];
            Start_Y=[563,1172,861,869,839,902,1573,2550,1512];
            Start_Z=[80,113,108,110,111,112,116,0,115];
        elseif strcmp(stack,'Section_03')
            Start_X=[2341,2094,1001,967,760,770];
            Start_Y=[286,790,1336,1807,1864,2001];
            Start_Z=[117,101,55,119,57,60];
        elseif strcmp(stack,'Section_04')
            Start_X=[2478,2487,3024,2417,2465,2443,3373,3513,2045,1452,2938,589,866,3262];
            Start_Y=[702,726,336,823,723,846,1037,1382,547,832,1533,1871,2066,1263];
            Start_Z=[23,84,16,88,16,19,15,20,5,22,25,26,30,92];
        elseif strcmp(stack,'Section_05')
            Start_X=[2710,3534,2725,1885,1152,336];
            Start_Y=[281,1073,1126,703,1174,1985];
            Start_Z=[6,10,13,5,12,13];
        end
        
        Start_X=(Start_X(start_point)/reduct_factor);
        Start_Y=(Start_Y(start_point)/reduct_factor);
        Start_Z=(Start_Z(start_point)/reduct_factor);      
       
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
        pth='C:\Armen\DIADEM\Data\Neuromuscular\Image Stacks\';
        if ~strcmp(stack,'mosaic')
            pth=[pth,stack,'\'];
            [Orig,sizeOrig,classOrig]=ImportStack(pth,data_format,relative_thr,reduct_type,reduct_factor,reduct_method);
            
            Start_X=200/reduct_factor; 
            Start_Y=1/reduct_factor; 
            Start_Z=20/reduct_factor;

        elseif strcmp(stack,'mosaic')
            load([pth,'\xyz_Aligned_Mosaic_3Median_0_005.mat']); 
            reduct_factor=1;
            data_format='Sparse';
%             Start_X=(SP(start_point,1)/reduct_factor);
%             Start_Y=(SP(start_point,2)/reduct_factor);
%             Start_Z=(SP(start_point,3)/reduct_factor);
        end
        
        S=1; % Center Surround Filter size
        thr_low=5000; % Threshold
        Region_size_thr=200; % Eliminate small regions
        
        Parameters.W=1; % Propagate the wave into all or one region: W=1,0
        Parameters.ShortTerminalBranchLength=15; % Remove terminal branches that are shorter than: ShortTerminalBranchLength
        Parameters.LongBranchLength=20; % Eliminate loops that are smaller than: LongBranchLength
        Parameters.ShortIntermediateBranchLength=20; % Remove intermediate branches that are shorter than: ShortIntermediateBranchLength
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
        Parameters.SmallTreeLength=30;
   
        
    case 'V1'
        pth='C:\Armen\DIADEM\Data\Visual Cortex\';
        if ~strcmp(stack,'mosaic')
            pth=[pth,stack,'\'];
            [Orig,sizeOrig,classOrig]=ImportStack(pth,data_format,relative_thr,reduct_type,reduct_factor,reduct_method);
            Orig=255-Orig;
            Start_X=1; 
            Start_Y=1; 
            Start_Z=1;

        elseif strcmp(stack,'mosaic')
            load([pth,'\xyz_Aligned_Mosaic_3Median_0_005.mat']); 
            reduct_factor=1;
            data_format='Sparse';
%             Start_X=(SP(start_point,1)/reduct_factor);
%             Start_Y=(SP(start_point,2)/reduct_factor);
%             Start_Z=(SP(start_point,3)/reduct_factor);
        end
        
        S=2; % Center Surround Filter size
        thr_low=30; % Threshold
        Region_size_thr=200; % Eliminate small regions
        
        Parameters.W=1; % Propagate the wave into all or one region: W=1,0
        Parameters.ShortTerminalBranchLength=15; % Remove terminal branches that are shorter than: ShortTerminalBranchLength
        Parameters.LongBranchLength=20; % Eliminate loops that are smaller than: LongBranchLength
        Parameters.ShortIntermediateBranchLength=20; % Remove intermediate branches that are shorter than: ShortIntermediateBranchLength
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

%Orig=Orig(200:end,200:end,:);
if strcmp(data_format,'Full')
%     T=0.001; % 0.01-0.001
%     alpha=0.5; % 0.5 % strength of the double-well potential
%     h=0.2; % threshold of the double-well potential
%     Im=Phase_Field(Orig,T,alpha,h);
    
    %Im = ReduceBackground(Orig,10);
    Im=CSF_cos(Orig,S); 
    Im=Threshold(Im,thr_low); 
    Im=EliminateSmallReg(Im,Region_size_thr);
    [AMlbl_merged,r_merged]=AutomatedTracing(Im,Start_X,Start_Y,Start_Z,Parameters);
elseif strcmp(data_format,'Sparse')
    Im=CenterSurroundFilter_Sparse(Orig,sizeOrig,S);
    Im=Threshold_Sparse(Im,thr_low);
    Im=EliminateSmallReg_Sparse(Im,sizeOrig,Region_size_thr);
    [AMlbl_merged,r_merged]=AutomatedTracing_Sparse(Im,sizeOrig,Start_X,Start_Y,Start_Z,Parameters);
end


%%%% Automated Tracing Parameters Structure %%%%
% 1. Propagate the wave into all or one region: W=1,0
% 2. Remove terminal branches that are shorter than: ShortTerminalBranchLength 
% 3. Eliminate loops that are smaller than: LongBranchLength
% 4. Remove intermediate branches that are shorter than: ShortIntermediateBranchLength
% 5. Trim branches near the xy faces of the stack and do not merge branches in these regions: trim
% 6. Number of points per um in the trace: pointsperum
% 7. Number of trace optimization steps: Nstep
% 8. Stiffness of the trace: alpha
% 9. Trace optimization step size: betta
% 10. Disconnect all branches at branch points: Disconnect=1,0
% 11. Consider merging branches that are closer than: MaxDistance
% 12. Number of steps taken to determine branch tip directions: stepsback
% 13. Minimal distance between branch tip clusters: Dist_thr
% 14. Merger cost function parameters: Cos2, d2, c2, m2, f2; Cos3, d3, c3, m3, f3; Cos4, d4, c4, m4, f4; g


% 3. Ploting the trace on a maximum projection image
figure
imshow(max(Im,[],3),[0 max(Im(:))]) 
hold on 
PlotAM(AMlbl_merged, r_merged)

%figure
% if strcmp(data_format,'Full')
%     imshow(max(Orig,[],3),[0 max(Orig(:))])
% elseif strcmp(data_format,'Sparse')
%     PlotMaxProjection_Sparse([],Orig,sizeOrig)
% end

% figure
PlotAM(AMlbl_merged, r_merged)
hold on
plot3(r_merged(1,2),r_merged(1,1),r_merged(1,3),'r*')
for i=1:1000
    i=mod(i,size(Orig,3));
    if i==0
        i=size(Orig,3);
    end
    i
    h=imshow(Orig(:,:,i),[0 max(Orig(:))]);
    pause
    delete(h)
    hold all
end


% 4. Exporting the trace
%[swc_trace_1,swc_trace_all] = AM2swc(AMlbl_merged,r_merged,reduct_type,reduct_factor); 
% Exporting reduced trace
%[swc_trace_1,swc_trace_all] = AM2swc(AMlbl_merged,[(r_merged(:,1:2)-1).*n,(r_merged(:,3)-1).*n]+1); 

% figure(2)
% imshow(max(Orig,[],3),[0 max(Orig(:))]);
% hold all; drawnow
% PlotSWCtrees(swc_trace_1)
%dlmwrite('C:\Armen\DIADEM\Data\Automated_Reconstructions\Neocortical_Layer_6_Axons\NC_01_1.swc', swc_trace_1, 'delimiter', ' ')
%dlmwrite('C:\Armen\DIADEM\Data\Automated_Reconstructions\Neocortical_Layer_6_Axons\NC_01_all.swc', swc_trace_all, 'delimiter', ' ')
%dlmwrite('C:\Armen\DIADEM\Data\Automated_Reconstructions\Olfactory_Projection_Fibers\OP_4.swc', swc_trace_1, 'delimiter', ' ')
%dlmwrite('C:\Armen\DIADEM\Data\Automated_Reconstructions\CA3\Section_02\9.swc', swc_trace_1, 'delimiter', ' ')
%dlmwrite('C:\Armen\DIADEM\Data\Automated_Reconstructions\Cerebellar Climbing Fibers\CF_2.swc', swc_trace_1, 'delimiter', ' ')




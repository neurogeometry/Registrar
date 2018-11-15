function Run_Dataset_Training(dataset,stack_num)

% L6: 01-06 (Qualifying) + 01-10 (Final) = 1-16
% Olfactory: OP_1-OP_3 + OP_4-OP_6 (Qualifying) + OP_1-OP_3 (Final) = 1-9
% Cerebellar: CF_1 + CF_2 (Qualifying) + CF_3 (Final) = 1-3
% CA3: Section_01 + Section_02-Section_05 (Qualifying) + Section 1-Section 4 (Final) = 1-9
% Neuromuscular: 001-152 (Qualifying) + 001-156 (Final) = 1-308
% Visual: Section 1-Section 4 (Final) = 1-53

data_format='Full'; 
relative_thr=0.1;
reduct_type='NA'; % 'NA','xy','xyz'
reduct_factor=1;
reduct_method='Median';
Start_X=1; Start_Y=1; Start_Z=1;

% 1. Import thesholded image stack
switch dataset
    case 'L6'
        pth1='C:\Armen\DIADEM\Qualifying Round Data\L6\Image Stacks\';
        pth2='C:\Armen\DIADEM\Final Round Data\L6\Image Stacks\';
        if stack_num>=1 && stack_num<=6
            stack_name=['0',num2str(stack_num)];
            pth=[pth1,stack_name,'\'];
        elseif stack_num>6 && stack_num<=15
            stack_name=['0',num2str(stack_num-6)];
            pth=[pth2,stack_name,'\'];
        elseif stack_num==16
            stack_name='10';
            pth=[pth2,stack_name,'\'];    
        else
            error('Stack numbers for this dataset range from 1 to 16')
        end
        [Orig,sizeOrig,classOrig]=ImportStack(pth,data_format,relative_thr,reduct_type,reduct_factor,reduct_method);
               
        S=1; % Center Surround Filter size 1,2
        thr_low=25; % Threshold 25,8
        Region_size_thr=200; % Eliminate small regions
        
        Parameters.W=1; % Propagate the wave into all or one region: W=1,0
        Parameters.ShortTerminalBranchLength=10; % Remove terminal branches that are shorter than: ShortTerminalBranchLength
        Parameters.LongBranchLength=20; % Eliminate loops that are smaller than: LongBranchLength
        Parameters.ShortIntermediateBranchLength=10; % Remove intermediate branches that are shorter than: ShortIntermediateBranchLength
        Parameters.trim=0; % Cut the trace near xy the faces of the stack
        Parameters.pointsperum=0.25; % Number of points per um in the trace: pointsperum
        Parameters.Nstep=20; % Number of trace optimization steps: Nstep
        Parameters.alpha=0.5; % Stiffness of the trace: alpha
        Parameters.betta=0.05; % Trace optimization step size: betta
        Parameters.sig_ips=1; % Intermediate-point sigma
        Parameters.sig_bps=2; % Branch-point sigma
        Parameters.sig_tps=2; % Terminal-point sigma
        Parameters.Disconnect=0; % Disconnect all branches at branch points: Disconnect=1,0
        Parameters.MaxDistance=17; % Consider merging branches that are closer than: MaxDistance
        Parameters.stepsback=40; % Number of steps taken to determine branch tip directions: stepsback
        Parameters.Dist_thr=17; % Minimal distance between branch tip clusters: Dist_thr
        Parameters.SmallTreeLength=30; % Eliminate small trees
        Parameters.R_min=1.0; % minimal neurite radius
        Parameters.R_step=0.5; % step in neurite radius
        Parameters.R_max=2; % maximal neurite radius
        Parameters.NBestScenarios=20; % # of best merger scenarios considered
        Parameters.FourWayMerger=0; % Allow 4-way branching or not: 1 or 0 
        
        % Initial cost function parameters:
        Parameters.Cos2=-1; Parameters.d2=0.1; Parameters.c2=1; Parameters.o2=-0.1; Parameters.m2=0; Parameters.f2=0;
        Parameters.Cos3=-0.5; Parameters.d3=0.1; Parameters.c3=1; Parameters.o3=-0.1; Parameters.m3=0; Parameters.f3=0;
        Parameters.Cos4=0; Parameters.d4=0.1; Parameters.c4=1; Parameters.o4=-0.1; Parameters.m4=0; Parameters.f4=0;
        Parameters.g=1.5;
        
        % Final cost function parameters:
%         Parameters.Ltotal=0.022284; Parameters.Itotal=-0.052874; Parameters.Imean=-0.00054294; 
%         Parameters.Tmax=0.64773; Parameters.CV=0.0032823; Parameters.Free_tips=0.22756; 
%         Parameters.Kmax=0.23124; Parameters.Kmean=0.68694;
        
        Parameters.Ltotal=0.0095186; Parameters.Itotal=-0.038776; Parameters.Imean=-0.080969; 
        Parameters.Tmax=0.025498; Parameters.CV=6.1029e-005; Parameters.Free_tips=0.049954; 
        Parameters.Kmax=0.55771; Parameters.Kmean=0.82321;
        
    case 'Olfactory'
        pth1='C:\Armen\DIADEM\Qualifying Round Data\Olfactory Projection Fibers\Image Stacks\Training Image Stack\';
        pth2='C:\Armen\DIADEM\Qualifying Round Data\Olfactory Projection Fibers\Image Stacks\Qualifier Image Stack\';
        pth3='C:\Armen\DIADEM\Final Round Data\Olfactory_Projection_Fibers\Image Stacks\';
        %Start_X=31.015; Start_Y=429.54; Start_Z=0;  
        %Start_X=0.72501; Start_Y=391.08; Start_Z=25;
        %Start_X=93.742; Start_Y=179; Start_Z=38;
        Start_X=128.2; Start_Y=504.37; Start_Z=0.303;
        %Start_X=185.7; Start_Y=264.02; Start_Z=33;
        %Start_X=15.074; Start_Y=412.01; Start_Z=10;
        if stack_num>=1 && stack_num<=3
            stack_name=['OP_',num2str(stack_num)];
            pth=[pth1,stack_name,'\'];
        elseif stack_num>3 && stack_num<=6
            stack_name=['OP_',num2str(stack_num)];
            pth=[pth2,stack_name,'\'];
        elseif stack_num>6 && stack_num<=9
            stack_name=['OP_',num2str(stack_num-6)];
            pth=[pth3,stack_name,'\'];    
        else
            error('Stack numbers for this dataset range from 1 to 9')
        end
        [Orig,sizeOrig,classOrig]=ImportStack(pth,data_format,relative_thr,reduct_type,reduct_factor,reduct_method);

        S=[1,1.5,2]; % Center Surround Filter size
        thr_low=7; % Threshold
        Region_size_thr=200; % Eliminate small regions
        
        Parameters.W=1; % Propagate the wave into all or one region: W=1,0
        Parameters.ShortTerminalBranchLength=8; % Remove terminal branches that are shorter than: ShortTerminalBranchLength
        Parameters.LongBranchLength=15; % Eliminate loops that are smaller than: LongBranchLength
        Parameters.ShortIntermediateBranchLength=6; % Remove intermediate branches that are shorter than: ShortIntermediateBranchLength
        Parameters.trim=0; % Cut the trace near the xy faces of the stack
        Parameters.pointsperum=0.25; % Number of points per um in the trace: pointsperum
        Parameters.Nstep=20; % Number of trace optimization steps: Nstep
        Parameters.alpha=0.5; % Stiffness of the trace: alpha
        Parameters.betta=0.05; % Trace optimization step size: betta
        Parameters.sig_ips=1; % Intermediate-point sigma
        Parameters.sig_bps=2; % Branch-point sigma
        Parameters.sig_tps=2; % Terminal-point sigma
        Parameters.Disconnect=0; % Disconnect all branches at branch points: Disconnect=1,0
        Parameters.MaxDistance=8; % Consider merging branches that are closer than: MaxDistance
        Parameters.stepsback=40; % Number of steps taken to determine branch tip directions: stepsback
        Parameters.Dist_thr=8; % Minimal distance between branch tip clusters: Dist_thr
        Parameters.SmallTreeLength=10; % Eliminate small trees
        
        % Initial cost function parameters:
        Parameters.Cos2=-1; Parameters.d2=0.1; Parameters.c2=2; Parameters.m2=-1; Parameters.f2=10;
        Parameters.Cos3=-0.5; Parameters.d3=0.1; Parameters.c3=0.5; Parameters.m3=-1; Parameters.f3=0;
        Parameters.Cos4=0; Parameters.d4=0.1; Parameters.c4=0; Parameters.m4=-3; Parameters.f4=0;
        Parameters.g=1;
        
        % Final cost function parameters:
        Parameters.Ltotal=0.0027153; Parameters.Itotal=-0.0057336; Parameters.Imean=-0;
        Parameters.Tmax=0.12354; Parameters.CV=0.00006; Parameters.Free_tips=0.024483;
        Parameters.Kmax=0.62279; Parameters.Kmean=0.77216;
        
    case 'Cerebellar'
        pth1='C:\Armen\DIADEM\Qualifying Round Data\Cerebellar Climbing Fibers\Image Stacks\Training Image Stack\';
        pth2='C:\Armen\DIADEM\Qualifying Round Data\Cerebellar Climbing Fibers\Image Stacks\Qualifier Image Stack\';
        pth3='C:\Armen\DIADEM\Final Round Data\Cerebellar Climbing Fibers\Image Stacks\';

        if stack_num==1
            stack_name='CF_1';
            pth=[pth1,stack_name,'\'];
            Start_X=4882; Start_Y=1797; Start_Z=19;
        elseif stack_num==2
            stack_name='CF_2';
            pth=[pth2,stack_name,'\'];
            Start_X=818; Start_Y=332; Start_Z=6;
        elseif stack_num==3
            stack_name='CF_3';
            pth=[pth3,stack_name,'\'];
            Start_X=1; Start_Y=1; Start_Z=1;
        else
            error('Stack numbers for this dataset range from 1 to 3')
        end
        
        if strcmp(reduct_type,'xy')
            Start_X=Start_X/reduct_factor;
            Start_Y=Start_Y/reduct_factor;
        elseif strcmp(reduct_type,'xyz')
            Start_X=Start_X/reduct_factor;
            Start_Y=Start_Y/reduct_factor;
            Start_Z=Start_Z/reduct_factor;
        end
        
        [Orig,sizeOrig,classOrig]=ImportStack(pth,data_format,relative_thr,reduct_type,reduct_factor,reduct_method);
        
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
        Parameters.SmallTreeLength=0;
        
        % Initial cost function parameters:
        Parameters.Cos2=-1; Parameters.d2=0.1; Parameters.c2=2; Parameters.m2=-1; Parameters.f2=0;
        Parameters.Cos3=-0.5; Parameters.d3=0.1; Parameters.c3=0.5; Parameters.m3=-1; Parameters.f3=0;
        Parameters.Cos4=0; Parameters.d4=100; Parameters.c4=1; Parameters.m4=0; Parameters.f4=0;
        Parameters.g=5;
                
        % Final cost function parameters:
        Parameters.Ltotal=0.0027153; Parameters.Itotal=-0.0057336; Parameters.Imean=-0;
        Parameters.Tmax=0.12354; Parameters.CV=0.00006; Parameters.Free_tips=0.024483;
        Parameters.Kmax=0.62279; Parameters.Kmean=0.77216;
        
    case 'CA3'
        pth1='C:\Armen\DIADEM\Qualifying Round Data\CA3\Training Image Stack\';
        pth2='C:\Armen\DIADEM\Qualifying Round Data\CA3\Qualifier Image Stacks\';
        pth3='C:\Armen\DIADEM\Final Round Data\CA3\Image Stacks\';
        if stack_num==1
            stack_name='Section_01';
            pth=[pth1,stack_name,'\'];
        elseif stack_num>1 && stack_num<=5
            stack_name=['Section_0',num2str(stack_num)];
            pth=[pth2,stack_name,'\'];
        elseif stack_num>5 && stack_num<=9
            stack_name=['Section ',num2str(stack_num-5)];
            pth=[pth3,stack_name,'\'];    
        else
            error('Stack numbers for this dataset range from 1 to 9')
        end
        [Orig,sizeOrig,classOrig]=ImportStack(pth,data_format,relative_thr,reduct_type,reduct_factor,reduct_method);
          
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
        pth1='C:\Armen\DIADEM\Qualifying Round Data\Neuromuscular\Image Stacks\';
        pth2='C:\Armen\DIADEM\Final Round Data\Neuromuscular\Image Stacks\';
        if stack_num>=1 && stack_num<=152
            stack_name=num2str(stack_num,'%03.0f');
            pth=[pth1,stack_name,'\'];
            [Orig,sizeOrig,classOrig]=ImportStack(pth,data_format,relative_thr,reduct_type,reduct_factor,reduct_method);
        elseif stack_num>152 && stack_num<=308
            stack_name=num2str(stack_num-152,'%03.0f');
            pth=[pth2,stack_name,'\'];  
            [Orig,sizeOrig,classOrig]=ImportStack(pth,data_format,relative_thr,reduct_type,reduct_factor,reduct_method);
            Orig=255-Orig; 
        else
            error('Stack numbers for this dataset range from 1 to 308')
        end
               
        S=2; % Center Surround Filter size
        thr_low=5000; % Threshold
        Region_size_thr=200; % Eliminate small regions
        
        Parameters.W=1; % Propagate the wave into all or one region: W=1,0
        Parameters.ShortTerminalBranchLength=15; % Remove terminal branches that are shorter than: ShortTerminalBranchLength
        Parameters.LongBranchLength=20; % Eliminate loops that are smaller than: LongBranchLength
        Parameters.ShortIntermediateBranchLength=20; % Remove intermediate branches that are shorter than: ShortIntermediateBranchLength
        Parameters.trim=0; % Cut the trace near the xy faces of the stack
        Parameters.pointsperum=0.25; % Number of points per um in the trace: pointsperum
        Parameters.Nstep=20; % Number of trace optimization steps: Nstep
        Parameters.alpha=0.5; % Stiffness of the trace: alpha
        Parameters.betta=0.1; % Trace optimization step size: betta
        Parameters.sig_ips=1; % Intermediate-point sigma
        Parameters.sig_bps=3; % Branch-point sigma
        Parameters.sig_tps=3; % Terminal-point sigma
        Parameters.Disconnect=0; % Disconnect all branches at branch points: Disconnect=1,0
        Parameters.MaxDistance=12; % Consider merging branches that are closer than: MaxDistance
        Parameters.stepsback=40; % Number of steps taken to determine branch tip directions: stepsback
        Parameters.Dist_thr=15; % Minimal distance between branch tip clusters: Dist_thr
        Parameters.SmallTreeLength=30; % Eliminate small trees
        
        % Initial cost function parameters:
        Parameters.Cos2=-1; Parameters.d2=0.1; Parameters.c2=2; Parameters.m2=-1; Parameters.f2=10;
        Parameters.Cos3=-0.5; Parameters.d3=0.2; Parameters.c3=2; Parameters.m3=-1; Parameters.f3=10;
        Parameters.Cos4=0; Parameters.d4=100; Parameters.c4=1; Parameters.m4=0; Parameters.f4=0;
        Parameters.g=2;
        
        % Final cost function parameters:
        Parameters.Ltotal=0.0024692; Parameters.Itotal=-0.024265; Parameters.Imean=-0.000025;
        Parameters.Tmax=0.020149; Parameters.CV=0.16013; Parameters.Free_tips=0.0093001;
        Parameters.Kmax=0.38378; Parameters.Kmean=0.90883;
        
    case 'Visual'
        pth1='C:\Armen\DIADEM\Final Round Data\Visual\Image Stacks\Section 1\';
        pth2='C:\Armen\DIADEM\Final Round Data\Visual\Image Stacks\Section 2\';
        pth3='C:\Armen\DIADEM\Final Round Data\Visual\Image Stacks\Section 3\';
        pth4='C:\Armen\DIADEM\Final Round Data\Visual\Image Stacks\Section 4\';
        if stack_num>=1 && stack_num<=6
            stack_names={'C01R1','C01R2','C01R3','C02R2','C02R3','C03R2'};
            pth=[pth1,stack_names{stack_num},'\']; 
        elseif stack_num>6 && stack_num<=16
            stack_names={'C01R2','C01R3','C02R3','C03R3','C04R1','C04R2','C05R1','C05R2','C06R1','C07R1'};
            pth=[pth2,stack_names{stack_num-6},'\'];  
        elseif stack_num>16 && stack_num<=38
            stack_names={'C01R3','C03R3','C04R1','C04R3','C05R1','C05R2','C05R3','C06R3','C07R1','C08R1','C08R2','C08R3','C09R1','C09R2','C10R1','C10R2','C11R1','C11R2','C12R1','C12R2','C13R1','C13R2'};
            pth=[pth3,stack_names{stack_num-16},'\'];      
        elseif stack_num>38 && stack_num<=53
            stack_names={'C01R2','C01R3','C03R3','C04R3','C10R1','C10R3','C11R1','C11R2','C11R3','C12R1','C12R2','C13R2','C14R2','C14R3','C15R2'};
            pth=[pth4,stack_names{stack_num-38},'\'];          
        else
            error('Stack numbers for this dataset range from 1 to 53')
        end
        [Orig,sizeOrig,classOrig]=ImportStack(pth,data_format,relative_thr,reduct_type,reduct_factor,reduct_method);
               
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

%Im=Threshold(Orig,0); %thr_low
%Im(Im==0)=120;
%Im=CSF_cos(Orig,S);
[Im, ~] = Multi_Scale_LoG(Orig,Parameters.R_min,Parameters.R_step,Parameters.R_max);
Im=Threshold(Im,thr_low); %
Im=Eliminate_Small_Regions(Im,Region_size_thr);
[AMlbl_initial,r_initial,R_initial,ClustersStr]=InitialTrace(Orig,Im,Start_X,Start_Y,Start_Z,Parameters);

figure
imshow(max(Orig,[],3),[0 max(Orig(:))]) 
hold on 
PlotAM(AMlbl_initial,r_initial)

%save('C:\Armen\Publications\Paper21 (Diadem)\Rebutal\Figures\CF Results\CF_1.mat','Orig','Im','AMlbl_merged','r_merged','R_merged','Parameters')
%[~,swc]=AM2swc(AMlbl_merged,r_merged,R_merged,reduct_type,reduct_factor); 
%dlmwrite('C:\Users\Armen Stepanyants\Desktop\Temp\Reconstructions\CF_1.swc',swc,'delimiter', ' ')
% for i=1:1000
%     i=mod(i,size(Orig,3));
%     if i==0
%         i=size(Orig,3);
%     end
%     i
%     h=imshow(Orig(:,:,i),[0 max(Orig(:))]);
%     pause
%     delete(h)
%     hold all
% end

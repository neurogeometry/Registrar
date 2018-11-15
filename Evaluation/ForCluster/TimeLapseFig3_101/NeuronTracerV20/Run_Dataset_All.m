function Run_Dataset_All(Data_Set,Stack_Num,Method)

reduction_x=1;
reduction_y=1;
reduction_z=1;
reduct_method='Max';

switch Data_Set
    case 'L6'
        pth1='C:\Armen\DIADEM\Qualifying Round Data\L6\Image Stacks\';
        pth2='C:\Armen\DIADEM\Final Round Data\L6\Image Stacks\';
        if Stack_Num>=1 && Stack_Num<=6
            stack_name=['0',num2str(Stack_Num)];
            pth=[pth1,stack_name,'\'];
        elseif Stack_Num>6 && Stack_Num<=15
            stack_name=['0',num2str(Stack_Num-6)];
            pth=[pth2,stack_name,'\'];
        elseif Stack_Num==16
            stack_name='10';
            pth=[pth2,stack_name,'\'];
        else
            error('Stack numbers for this dataset range from 1 to 16')
        end
                
        R_min=2; % minimal neurite radius
        R_step=1; % step in neurite radius
        R_max=3; % maximal neurite radius
        thr_low=12; % Threshold 25,8
        Region_size_thr=100; % Eliminate small regions
        N_steps=2000;
        Max_Known_Dist=15;
        Contrast=0.2;
        unisotropy=[1,1,1/5];
        
        parameter_file_path='C:\Armen\DIADEM\Neuron Tracer V19\Parameter Files\Parameters_L6.txt';
        TrainingHistory_file_path='C:\Armen\DIADEM\Neuron Tracer V19\Parameter Files\TrainingHistory_L6.mat';
        
    case 'Olfactory'
        pth1='C:\Armen\DIADEM\Qualifying Round Data\Olfactory Projection Fibers\Image Stacks\Training Image Stack\';
        pth2='C:\Armen\DIADEM\Qualifying Round Data\Olfactory Projection Fibers\Image Stacks\Qualifier Image Stack\';
        pth3='C:\Armen\DIADEM\Final Round Data\Olfactory_Projection_Fibers\Image Stacks\';

        if Stack_Num>=1 && Stack_Num<=3
            stack_name=['OP_',num2str(Stack_Num)];
            pth=[pth1,stack_name,'\'];
        elseif Stack_Num>3 && Stack_Num<=6
            stack_name=['OP_',num2str(Stack_Num)];
            pth=[pth2,stack_name,'\'];
        elseif Stack_Num>6 && Stack_Num<=9
            stack_name=['OP_',num2str(Stack_Num-6)];
            pth=[pth3,stack_name,'\'];
        else
            error('Stack numbers for this dataset range from 1 to 9')
        end
        
        R_min=2; % minimal neurite radius
        R_step=2; % step in neurite radius
        R_max=4; % maximal neurite radius
        thr_low=7; % Threshold
        Region_size_thr=200; % Eliminate small regions
        N_steps=2000;
        Max_Known_Dist=10;
        Contrast=0.2;
        unisotropy=[1,1,1/5];
        
        parameter_file_path='C:\Armen\DIADEM\Neuron Tracer V19\Parameter Files\Parameters_OP.txt';
        TrainingHistory_file_path='C:\Armen\DIADEM\Neuron Tracer V19\Parameter Files\TrainingHistory_OP.mat';
        
    case 'Cerebellar'
        pth1='C:\Armen\DIADEM\Qualifying Round Data\Cerebellar Climbing Fibers\Image Stacks\Training Image Stack\';
        pth2='C:\Armen\DIADEM\Qualifying Round Data\Cerebellar Climbing Fibers\Image Stacks\Qualifier Image Stack\';
        pth3='C:\Armen\DIADEM\Final Round Data\Cerebellar Climbing Fibers\Image Stacks\';
        
        if Stack_Num==1
            stack_name='CF_1';
            pth=[pth1,stack_name,'\'];
        elseif Stack_Num==2
            stack_name='CF_2';
            pth=[pth2,stack_name,'\'];
        elseif Stack_Num==3
            stack_name='CF_3';
            pth=[pth3,stack_name,'\'];
        else
            error('Stack numbers for this dataset range from 1 to 3')
        end
        
        R_min=2; % minimal neurite radius
        R_step=2; % step in neurite radius
        R_max=2; % maximal neurite radius
        thr_low=20; % Threshold 25,8
        Region_size_thr=200; % Eliminate small regions
        N_steps=2000;
        Max_Known_Dist=15;
        Contrast=0.2;
        unisotropy=[1,1,3];
        
    case 'CA3'
        pth1='C:\Armen\DIADEM\Qualifying Round Data\CA3\Training Image Stack\';
        pth2='C:\Armen\DIADEM\Qualifying Round Data\CA3\Qualifier Image Stacks\';
        pth3='C:\Armen\DIADEM\Final Round Data\CA3\Image Stacks\';
        if Stack_Num==1
            stack_name='Section_01';
            pth=[pth1,stack_name,'\'];
        elseif Stack_Num>1 && Stack_Num<=5
            stack_name=['Section_0',num2str(Stack_Num)];
            pth=[pth2,stack_name,'\'];
        elseif Stack_Num>5 && Stack_Num<=9
            stack_name=['Section ',num2str(Stack_Num-5)];
            pth=[pth3,stack_name,'\'];
        else
            error('Stack numbers for this dataset range from 1 to 9')
        end
        
        R_min=0; % minimal neurite radius
        R_step=0; % step in neurite radius
        R_max=0; % maximal neurite radius
        thr_low=50; % Threshold
        Region_size_thr=200; % Eliminate small regions
        N_steps=2000;
        Max_Known_Dist=15;
        Contrast=0.2;
        unisotropy=[1,1,1/5];
        
        parameter_file_path='C:\Armen\DIADEM\Neuron Tracer V19\Parameter Files\Parameters_OP.txt';
        TrainingHistory_file_path='C:\Armen\DIADEM\Neuron Tracer V19\Parameter Files\TrainingHistory_OP.mat';
        
        
    case 'Neuromuscular'
        pth1='C:\Armen\DIADEM\Qualifying Round Data\Neuromuscular\Image Stacks\';
        pth2='C:\Armen\DIADEM\Final Round Data\Neuromuscular\Image Stacks\';
        if Stack_Num>=1 && Stack_Num<=152
            stack_name=num2str(Stack_Num,'%03.0f');
            pth=[pth1,stack_name,'\'];
        elseif Stack_Num>152 && Stack_Num<=308
            stack_name=num2str(Stack_Num-152,'%03.0f');
            pth=[pth2,stack_name,'\'];
        else
            error('Stack numbers for this dataset range from 1 to 308')
        end
        
        R_min=4; % minimal neurite radius
        R_step=2; % step in neurite radius
        R_max=10; % maximal neurite radius
        thr_low=5000; % Threshold
        Region_size_thr=200; % Eliminate small regions
        N_steps=2000;
        Max_Known_Dist=30;
        Contrast=0.2;
        unisotropy=[1,1,1];
        
        parameter_file_path='C:\Armen\DIADEM\Neuron Tracer V19\Parameter Files\Parameters_OP.txt';
        TrainingHistory_file_path='C:\Armen\DIADEM\Neuron Tracer V19\Parameter Files\TrainingHistory_OP.mat';
        
    case 'Visual'
        pth1='C:\Armen\DIADEM\Final Round Data\Visual Cortex\Image Stacks\Section 1\';
        pth2='C:\Armen\DIADEM\Final Round Data\Visual Cortex\Image Stacks\Section 2\';
        pth3='C:\Armen\DIADEM\Final Round Data\Visual Cortex\Image Stacks\Section 3\';
        pth4='C:\Armen\DIADEM\Final Round Data\Visual Cortex\Image Stacks\Section 4\';
        if Stack_Num>=1 && Stack_Num<=6
            stack_names={'C01R1','C01R2','C01R3','C02R2','C02R3','C03R2'};
            pth=[pth1,stack_names{Stack_Num},'\'];
        elseif Stack_Num>6 && Stack_Num<=16
            stack_names={'C01R2','C01R3','C02R3','C03R3','C04R1','C04R2','C05R1','C05R2','C06R1','C07R1'};
            pth=[pth2,stack_names{Stack_Num-6},'\'];
        elseif Stack_Num>16 && Stack_Num<=38
            stack_names={'C01R3','C03R3','C04R1','C04R3','C05R1','C05R2','C05R3','C06R3','C07R1','C08R1','C08R2','C08R3','C09R1','C09R2','C10R1','C10R2','C11R1','C11R2','C12R1','C12R2','C13R1','C13R2'};
            pth=[pth3,stack_names{Stack_Num-16},'\'];
        elseif Stack_Num>38 && Stack_Num<=53
            stack_names={'C01R2','C01R3','C03R3','C04R3','C10R1','C10R3','C11R1','C11R2','C11R3','C12R1','C12R2','C13R2','C14R2','C14R3','C15R2'};
            pth=[pth4,stack_names{Stack_Num-38},'\'];
        else
            error('Stack numbers for this dataset range from 1 to 53')
        end
        
        R_min=0; % minimal neurite radius
        R_step=0; % step in neurite radius
        R_max=0; % maximal neurite radius
        thr_low=30; % Threshold
        Region_size_thr=200; % Eliminate small regions
        N_steps=2000;
        Max_Known_Dist=15;
        Contrast=0.2;
        unisotropy=[1,1,3];
        
    case 'Kathy'
        pth='C:\Armen\Kathy\63x hi res scans\090318_R188_138\ch00\';
        
        R_min=0; % minimal neurite radius
        R_step=0; % step in neurite radius
        R_max=0; % maximal neurite radius
        thr_low=25; % Threshold 25,8
        Region_size_thr=200; % Eliminate small regions
        N_steps=2000;
        Max_Known_Dist=15;
        Contrast=0.2;
        unisotropy=[1,1,3];
        
    case 'Test'
        %pth='C:\Users\Armen Stepanyants\Desktop\Kathy_139_B\';
        %pth='C:\Armen\Uchida\calvGATSNC3\';
        %pth='C:\Users\Armen Stepanyants\Desktop\Green2\';
        %pth='C:\Users\Armen Stepanyants\Desktop\Images\';
        pth='C:\Users\Armen Stepanyants\Desktop\Images\';
                
        R_min=2; % minimal neurite radius
        R_step=1; % step in neurite radius
        R_max=4; % maximal neurite radius
        thr_low=7; % Threshold
        Region_size_thr=200; % Eliminate small regions
        N_steps=2000;
        Max_Known_Dist=15;
        Contrast=0.2;
        unisotropy=[1,1,3];
        
        parameter_file_path='C:\Armen\DIADEM\Neuron Tracer V19\Parameter Files\Parameters_OP.txt';
        TrainingHistory_file_path='C:\Armen\DIADEM\Neuron Tracer V19\Parameter Files\TrainingHistory_OP.mat';
        
    case 'AFM'
        pth='C:\Armen\Nathan\AFM_1\';
        
        R_min=3; % minimal neurite radius
        R_step=1; % step in neurite radius
        R_max=3; % maximal neurite radius
        thr_low=20; % Threshold
        Region_size_thr=200; % Eliminate small regions
        N_steps=2000;
        Max_Known_Dist=15;
        Contrast=0.2;
        unisotropy=[1,1,3];
        
    case 'Karel'
        pth='C:\Armen\Karel\00129\';
        
        R_min=3; % minimal neurite radius
        R_step=1; % step in neurite radius
        R_max=3; % maximal neurite radius
        thr_low=20; % Threshold
        Region_size_thr=200; % Eliminate small regions
        N_steps=2000;
        Max_Known_Dist=15;
        Contrast=0.2;
        unisotropy=[1,1,3];
        
    case 'Steve'
        pth='C:\Armen\SteveAmato\Neuron 2R\';
        
        R_min=3; % minimal neurite radius
        R_step=1; % step in neurite radius
        R_max=3; % maximal neurite radius
        thr_low=20; % Threshold
        Region_size_thr=200; % Eliminate small regions
        N_steps=2000;
        Max_Known_Dist=15;
        Contrast=0.2;
        unisotropy=[1,1,3];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Load the stack
temp=dir([pth,'\*.tif']);
Names={temp.name};
for i=1:length(Names)
    Names{i}=(Names{i}(1:find(Names{i}=='.')-1));
end
[~,ind]=sort(str2double(Names));
Names=Names(ind);
for i=1:length(Names)
    Names{i}=[Names{i},'.tif'];
end

[Orig,~,~]=ImportStackJ(pth,Names,reduction_x,reduction_y,reduction_z,reduct_method);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Load parameters
Parameters=Read_Parameter_File(parameter_file_path);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Automated tracing

if strcmp(Method,'FM')
    Im=Orig;
    SVr=Find_Seeds(Orig,25,50,R_min,R_step,R_max,[]);
    [AMlbl,r,R]=InitialTrace_FM(Orig,Im,SVr,N_steps,Max_Known_Dist,Contrast,R_min,R_step,R_max,unisotropy,Parameters);
    
elseif strcmp(Method,'VC')
    Im = Multi_Scale_LoG(Orig,R_min,R_step,R_max);
    Im=Threshold_Background(Im,thr_low);
    Im=Eliminate_Small_Regions(Im,Region_size_thr);
    [AMlbl,r,R]=AutomatedTracing(Orig,Im,Parameters,TrainingHistory_file_path);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure
imshow(max(Orig,[],3),[0 max(Orig(:))])
hold on
PlotAM(AMlbl, r)

figure
imshow(max(Orig,[],3),[0 max(Orig(:))])
hold on
PlotR(AMlbl, r, R)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Save the results in the NCTracer format

%filepath=['C:\Armen\Publications\Paper24 (active learning)\Reconstructions_Auto_FM2\',Data_Set,'_',num2str(Stack_Num),'.mat'];
%Original=Orig; IM=Original; AM=AMlbl;
%save(filepath,'Original','IM','AM','r','R','reduction_xy','reduction_z')

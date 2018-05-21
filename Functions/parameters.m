%% This file store all available parameters
                            
                                %% Feature Extraction Parameters

    % Parameter to set the filter value
params.FE.filterValue = 1;          % For Smooth3 set to 1, 3D Gaussian set to 2, Imboxfilter3 set to 3, MultiScale LOG set to 4

    % Parameter to control the smooth3 filtering size
params.FE.smooth3BoxSize = 3;       %default:   params.FE.smooth3BoxSize = 3

    % Parameter to control the Gaussian filtering size
params.FE.GaussianSize = 3;         %default:   params.FE.GaussianSize = 3

    % Parameter to control the imbox filtering size
params.FE.IMboxSize = [5 5 3];      %default:   params.FE.IMboxSize = [5 5 3]

    % Parameter to control the Multi Scale LOG filtering size
params.FE.MLOG1 = 3;               %default:   params.FE.MLOG1 = 3
params.FE.MLOG2 = 3;               %default:   params.FE.MLOG2 = 3
params.FE.MLOG3 = 3;               %default:   params.FE.MLOG3 = 3

    %  Size of the window to exgtract information around features 
params.FE.windowsize = [4 4 4];

    % Expected number of features per stack
params.FE.MaxN_Features=1000;
params.FE.Mesh = [15 15 5]; %Features distance / Higher for less Features

    % Expected [x,y,z] movement of voxles during registration relative to stack size
    % Use 1 if unknown 
params.FE.Expected_Missalignment = [0.1,0.1,0.1];                                
                                %% Feature Matching Parameters
                 
    % Parameter controling the minimum matched points to run Ransac in Stitching_3D_Func.m 
paramsFMRansacMin =15;             %default:   params.FM.RansacLimit=15

    % Parameter controling the maximum distance of features to be consider in Hungarian Algorithm in Stitching_3D_Func.m 
paramsFMDT = 60;                   %default:	params.FM.DT = 60;       


    % Parameter to selec the similarity metric for feature matching
% available metrics includes:
% ZMAD = ZeroMeanAbsoluteDifferences
% MSE: MeanSquareError
% AD = AverageDifferences
% NAE = NormalizedAbsoluteError
% NCC= NormalizedCrossCorrelation
% ZNCC= Zero mean NormalizedCrossCorrelation
% SAD = Sum of Absolute Differences
% SAD =  Sum of Squared Differences
paramsFMMetric = 'ZMAD';        %default:	paramsFMMetric = 'ZMAD';

    % Parameter to control the maximum iteration of Ransac
paramsFMmaxiter=1000;           %default:	paramsFMmaxiter = 1000;

% Parameter to control the debug mode
% paramsFMdebug=0;

    % Use Hessian-Eigen
paramsFMEigen=0;                %default:	paramsFMEigen = 0;

% Minimum required features to do feature matching
% paramsFMminReqFeatures = 15;

params.FM.FeaturesFolderName='/features/';
params.FM.MatchedLocationsFile_Translation='/MatchedPoints_Translation.mat';
params.FM.MatchedLocationsFile_Rigid = '/MatchedPoints_Rigid.mat';
params.FM.MatchedLocationsFile_Affine = '/MatchedPoints_Affine.mat';
params.FM.MatchedLocationsFile_Non_Rigid = '/MatchedPoints_Non-Rigid.mat';

params.FM.MatchedLocationsFile_Translation_Global='/MatchedPoints_Translation_Global.mat';
params.FM.MatchedLocationsFile_Rigid_Global = '/MatchedPoints_Rigid_Global.mat';
params.FM.MatchedLocationsFile_Affine_Global = '/MatchedPoints_Affine_Global.mat';
params.FM.MatchedLocationsFile_Non_Rigid_Global = '/MatchedPoints_Non-Rigid_Global.mat';

% params.FM.AffineTransformFile='/Affine.mat';
% params.FM.TranslationTransformFile='/Translation.mat';
params.FM.StackPositions_Registered1='/StackPositions_Registered.mat';


                                %% Global Optimization Parameters
%Parameters related to generation of profiles used in gui_optimization.m and profilefilters.m
%regularization of individual shifts
params.GT.lamda=10^-6;
params.GA.lamda=10^-6;

    %regularization of the overall shift 
% params.GO.eps = 0.0095;
params.GT.eta = 1;
params.GA.eta = 1;

    % regularization of the overall deformation
params.GA.nu = 1;
    % regularization of individual deformations 
params.GA.mu=10^-1;

    % Registratio Results save file GT = Global Translation, GA = Global Affine
params.G.StackPositions_Registered = '\StackPositions_Registered';
params.G.StackPositions_RegisteredCSV = '\StackPositions_Registered.csv';

params.GT.StackPositions_Registered = '\StackPositions_Registered_Translation.mat';
params.GA.StackPositions_Registered = '\StackPositions_Registered_Affine.mat';
params.GR.StackPositions_Registered = '\StackPositions_Registered_Rigid.mat';
params.GN.StackPositions_Registered = '\StackPositions_Registered_nonrigid.mat';
params.GT.StackPositions_RegisteredCSV = '\StackPositions_Registered_Translation.csv';
params.GA.StackPositions_RegisteredCSV = '\StackPositions_Registered_Affine.csv';
params.GR.StackPositions_RegisteredCSV = '\StackPositions_Registered_Rigid.csv';
params.GN.StackPositions_RegisteredCSV = '\StackPositions_Registered_nonrigid.csv';

params.GT.Transformation = '\T_Translation.mat';
params.GR.Transformation = '\T_Rigid.mat';
params.GA.Transformation = '\T_Affine.mat';


                         %% Blending for Preview Parameters
%Parameters to control the size of the visualization tile size
params.BP.extSize = [300 300 0];

params.BP.saveImages = 0;
% Set to 1 if your stack has grey pad to remove in visualization
params.BP.remove_pad = 0;

                         %% Resampling Parameters
params.RE.savefolder = '\Tiles\';
paramsREuseHDF5 = 0;
paramsBigTileSize = [1024 1024 256]; % powers of 2 only
paramsFinalTileSize = [128 128 64]; % powers of 2 only
                         
                         
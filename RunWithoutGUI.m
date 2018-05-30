close all;
clear all;
clc

% Input CSV file
StackList_csv_pth = '../../MicroscopeFiles\Neocortical2_StackList.csv';
% Transformation Type, 1: Translation, 2: Rigid, 3: Affine, 4:Non-Rigid
TransformationValue = 1;
% Run in 1: Sequential or 2: Parallel
Seq_Par = 1;
% Number of Workes 
Par_workers = 6;

registeration (StackList_csv_pth,TransformationValue,Seq_Par,Par_workers);
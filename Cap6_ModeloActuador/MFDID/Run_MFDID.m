% Curve fitting using Dr. Kim's program (MFDID)
% coded by K. S. Park at UIUC

clear all; close all; clc;

fprintf('Run MFDID ...\n');
load pexdata        % load pseudo experimental data

TFw=[TFd1w_n TFa1w_n TFa2w_n TFa3w_n];      % transfer function from excitation (SIMO system)
TFu=[TFd1u_n TFa1u_n TFa2u_n TFa3u_n];      % transfer function from control input (SIMO system)
TF=[TFw TFu];                               % transfer function from excitation and control inputs (MIMO system)
save tmp TFw TFu TF freq
clear all; load tmp
mfdid                                       % run MFDID program
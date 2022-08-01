% Construct state-space model from curve fitting data
% coded by K. S. Park at UIUC

clear all; close all; clc;

fprintf('Get state-space model from results of MFDID ...\n');
% For transfer function from earthquake excitation (SIMO case)
load TFw_esti
[z,p,k]=zpkdata(sys_mfdid);
[Aew,Bew,Cew,Dew]=ssdata(sys_mfdid);
% check the stability
if (max(real(eig(Aew))) > 0)
    disp('unstable Ae')
end

% For transfer function from control input (SIMO case)
load TFu_esti
[z,p,k]=zpkdata(sys_mfdid);
[Aeu,Beu,Ceu,Deu]=ssdata(sys_mfdid);
% check the stability
if (max(real(eig(Aeu))) > 0)
    disp('unstable Ae')
end

% For transfer function from excitation and control inputs (MIMO case)
load TF_esti
[z,p,k]=zpkdata(sys_mfdid);
[Ae,Be,Ce,De]=ssdata(sys_mfdid);
% check the stability
if (max(real(eig(Ae))) > 0)
    disp('unstable Ae')
end
save ABCD_ESTI Aew Bew Cew Dew Aeu Beu Ceu Deu Ae Be Ce De
fprintf('Get state-space model from results of MFDID ...Done\n');



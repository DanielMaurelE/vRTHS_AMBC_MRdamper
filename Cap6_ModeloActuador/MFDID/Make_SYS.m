% Make state-space form of analytical model
% coded by K. S. Park at UIUC

clear all; close all; clc;
fprintf('Constructing analytical model for pseudo experimental data for MFDID ...\n');
% SYSTEM PARAMETERS
% L. L. Chung, R. C. Lin, T. T. Soong, A. M. Reinhorn,
% "Experimental study of active control for MDOF seismic structures",
% Journal of Engineering Mechanics, ASCE, 115(8), pp. 1609 ~ 1627 (1998).
mm=[5.6 0 0;0 5.6 0;0 0 5.6];                                      % lb-sec^2/in
kk=[15649 -9370 2107;-9370 17250 -9274;2107 -9274 7612];           % lb/in
cc=[2.185 -0.327 0.352;-0.327 2.608 -0.015;0.352 -0.015 2.497];    % lb-sec/in
ndof=size(mm,1);                                                   % number of degree of freedom

% ACTUATOR PARAMETERS
% S. J. Dyke, B. F. Spencer, Jr., P. Quast, M. K. Sain,
% "Role of control-structure interaction in protective system design",
% Journal of Engineering Mechanics, ASCE, 121(2), pp. 322-338 (1995)
% and Lecture note of Prof. B. F. Spencer, Jr.
% df=ca1*(u-x1)-ca2*dx1-ca3*f
ca1=2.0833e6;   % 2*beta*kq*A*gamma/V (lb/in-sec)
ca2=1.25e5;     % 2*beta*A^2/V (lb/in)
ca3=66.67;      % 2*beta*kc/V (1/sec)

% STATE-SPACE FORM REALIZATION
% w/ actuator
Bs=[1 0 0]';            % control force vector
% Gs=[0 0 1]';          % excitation vector for hammer test
Gs=-mm*ones(ndof,1);   % excitation vector for earthquake loading

A1=[zeros(ndof,ndof) eye(ndof) zeros(3,1); -inv(mm)*kk -inv(mm)*cc inv(mm)*Bs; ...
    -ca1 0 0 -ca2 0 0 -ca3];
Bu=[zeros(2*ndof,1);ca1];
Bg=[zeros(ndof,1);inv(mm)*Gs;0];
B1=[Bg Bu];
% for all state (x1 x2 x3 dx1 dx2 dx3 ddxa1 ddxa2 ddxa3 f)
C1=[eye(3) zeros(3,3) zeros(3,1);zeros(3,3) eye(3) zeros(3,1);-inv(mm)*kk -inv(mm)*cc inv(mm)*Bs;zeros(1,6) 1];
D1=zeros(size(C1,1),size(B1,2));

% w/o actuator
A2=A1(1:6,1:6);
B2=Bg(1:6,:);
C2=C1(1:9,1:6);
D2=D1(1:9,1);

save ABCD A1 B1 C1 D1 A2 B2 C2 D2
fprintf('Constructing analytical model for pseudo experimental data for MFDID ...Done !\n');
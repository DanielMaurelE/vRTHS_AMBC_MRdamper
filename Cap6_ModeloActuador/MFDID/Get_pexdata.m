% Get pseudo experimental data from analytical model
% coded by K. S. Park at UIUC

clear all; close all; clc;

fprintf('Get transfer function data from analytical model using frequency response function ...\n');
load ABCD                           % load system matrix
% TRANSFER FUNCTION
% set parameters
nfft=2048;                          % number of fft
fc=50;                              % cutoff frequency (Hz)
dt=1/(fc*2.56);                     % sampling time (sec)
fs=1/dt;                            % sampling frequency (Hz)
freq=[0:nfft]/nfft*(fc*1.28);       % frequency vector (Hz)
omega=[0:nfft]/nfft*(fc*1.28)*2*pi; % frequency vector (rad/sec)

% frequency response function for excitation input
% magnitude
TFd1w=freqresp(A1,B1,C1(1,:),D1(1,:),1,sqrt(-1)*omega);     % displacement of first story
TFa1w=freqresp(A1,B1,C1(7,:),D1(7,:),1,sqrt(-1)*omega);     % acceleration of first story
TFa2w=freqresp(A1,B1,C1(8,:),D1(8,:),1,sqrt(-1)*omega);     % acceleration of second story
TFa3w=freqresp(A1,B1,C1(9,:),D1(9,:),1,sqrt(-1)*omega);     % acceleration of third story
TFfw=freqresp(A1,B1,C1(10,:),D1(10,:),1,sqrt(-1)*omega);    % control force
% phase angle - need for comparison later
PAd1w=mod(angle(TFd1w)*180/pi+180,360)-180;
PAa1w=mod(angle(TFa1w)*180/pi+180,360)-180;
PAa2w=mod(angle(TFa2w)*180/pi+180,360)-180;
PAa3w=mod(angle(TFa3w)*180/pi+180,360)-180;
PAfw=mod(angle(TFfw)*180/pi+180,360)-180;
% add noise to emulate the experiment
N=size(TFd1w,1);
noise_percent=15;   % noise percent
tf_noise=(rand(N,1)-0.5)+sqrt(-1)*(rand(N,1)-0.5); 
% noised transfer function (magnitud and phase) - only magnitude will be used in MFDID 
TFd1w_n=TFd1w+tf_noise.*(noise_percent/100*abs(TFd1w)./abs(tf_noise));
TFa1w_n=TFa1w+tf_noise.*(noise_percent/100*abs(TFa1w)./abs(tf_noise));
TFa2w_n=TFa2w+tf_noise.*(noise_percent/100*abs(TFa2w)./abs(tf_noise));
TFa3w_n=TFa3w+tf_noise.*(noise_percent/100*abs(TFa3w)./abs(tf_noise));
TFfw_n=TFfw+tf_noise.*(noise_percent/100*abs(TFfw)./abs(tf_noise));
PAd1w_n=mod(angle(TFd1w_n)*180/pi+180,360)-180;
PAa1w_n=mod(angle(TFa1w_n)*180/pi+180,360)-180;
PAa2w_n=mod(angle(TFa2w_n)*180/pi+180,360)-180;
PAa3w_n=mod(angle(TFa3w_n)*180/pi+180,360)-180;
PAfw_n=mod(angle(TFfw_n)*180/pi+180,360)-180;

% plot the transfer functions
figure
subplot(221),plot(freq,db(abs(TFd1w)),freq,db(abs(TFd1w_n))),grid,xlabel('freq(Hz)'),ylabel('TFd1w (dB)'),legend('w/o noise','w/ noise')
subplot(223),plot(freq,PAd1w,freq,PAd1w_n),grid,xlabel('freq(Hz)'),ylabel('PAd1w (deg)')
subplot(222),plot(freq,db(abs(TFfw)),freq,db(abs(TFfw_n))),grid,xlabel('freq(Hz)'),ylabel('TFfw (dB)'),legend('w/o noise','w/ noise')
subplot(224),plot(freq,PAfw,freq,PAfw_n),grid,xlabel('freq(Hz)'),ylabel('PAfw (deg)')
figure
subplot(221),plot(freq,db(abs(TFa1w)),freq,db(abs(TFa1w_n))),grid,xlabel('freq(Hz)'),ylabel('TFa1w (dB)'),legend('w/o noise','w/ noise')
subplot(223),plot(freq,PAa1w,freq,PAa1w_n),grid,xlabel('freq(Hz)'),ylabel('PAa1w (deg)')
subplot(222),plot(freq,db(abs(TFa2w)),freq,db(abs(TFa2w_n))),grid,xlabel('freq(Hz)'),ylabel('TFa2w (dB)'),legend('w/o noise','w/ noise')
subplot(224),plot(freq,PAa2w,freq,PAa2w_n),grid,xlabel('freq(Hz)'),ylabel('PAa2w (deg)')
figure
subplot(221),plot(freq,db(abs(TFa3w)),freq,db(abs(TFa3w_n))),grid,xlabel('freq(Hz)'),ylabel('TFa3h (dB)'),legend('w/o noise','w/ noise')
subplot(223),plot(freq,PAa3w,freq,PAa3w_n),grid,xlabel('freq(Hz)'),ylabel('PAa3h (deg)')

% frequency response function for control input
% magnitude
TFd1u=freqresp(A1,B1,C1(1,:),D1(1,:),2,sqrt(-1)*omega);     % displacement of first story
TFa1u=freqresp(A1,B1,C1(7,:),D1(7,:),2,sqrt(-1)*omega);     % acceleration of first story
TFa2u=freqresp(A1,B1,C1(8,:),D1(8,:),2,sqrt(-1)*omega);     % acceleration of second story
TFa3u=freqresp(A1,B1,C1(9,:),D1(9,:),2,sqrt(-1)*omega);     % acceleration of third story
TFfu=freqresp(A1,B1,C1(10,:),D1(10,:),2,sqrt(-1)*omega);    % control force
% phase angle
PAd1u=mod(angle(TFd1u)*180/pi+180,360)-180;
PAa1u=mod(angle(TFa1u)*180/pi+180,360)-180;
PAa2u=mod(angle(TFa2u)*180/pi+180,360)-180;
PAa3u=mod(angle(TFa3u)*180/pi+180,360)-180;
PAfu=mod(angle(TFfu)*180/pi+180,360)-180;

% add noise to emulate the experiment
N=size(TFd1u,1);
tf_noise=(rand(N,1)-0.5)+sqrt(-1)*(rand(N,1)-0.5); 
% noised transfer function (magnitud and phase) - only magnitude will be used in MFDID 
TFd1u_n=TFd1u+tf_noise.*(noise_percent/100*abs(TFd1u)./abs(tf_noise));
TFa1u_n=TFa1u+tf_noise.*(noise_percent/100*abs(TFa1u)./abs(tf_noise));
TFa2u_n=TFa2u+tf_noise.*(noise_percent/100*abs(TFa2u)./abs(tf_noise));
TFa3u_n=TFa3u+tf_noise.*(noise_percent/100*abs(TFa3u)./abs(tf_noise));
TFfu_n=TFfu+tf_noise.*(noise_percent/100*abs(TFfu)./abs(tf_noise));
PAd1u_n=mod(angle(TFd1u_n)*180/pi+180,360)-180;
PAa1u_n=mod(angle(TFa1u_n)*180/pi+180,360)-180;
PAa2u_n=mod(angle(TFa2u_n)*180/pi+180,360)-180;
PAa3u_n=mod(angle(TFa3u_n)*180/pi+180,360)-180;
PAfu_n=mod(angle(TFfu_n)*180/pi+180,360)-180;

% plot the transfer functions
figure
subplot(221),plot(freq,db(abs(TFd1u)),freq,db(abs(TFd1u_n))),grid,xlabel('freq(Hz)'),ylabel('TFd1u (dB)'),legend('w/o noise','w/ noise')
subplot(223),plot(freq,PAd1u,freq,PAd1u_n),grid,xlabel('freq(Hz)'),ylabel('PAd1h (deg)')
subplot(222),plot(freq,db(abs(TFfu)),freq,db(abs(TFfu_n))),grid,xlabel('freq(Hz)'),ylabel('TFfu (dB)'),legend('w/o noise','w/ noise')
subplot(224),plot(freq,PAfu,freq,PAfu_n),grid,xlabel('freq(Hz)'),ylabel('PAfu (deg)')
figure
subplot(221),plot(freq,db(abs(TFa1u)),freq,db(abs(TFa1u_n))),grid,xlabel('freq(Hz)'),ylabel('TFa1u (dB)'),legend('w/o noise','w/ noise')
subplot(223),plot(freq,PAa1u,freq,PAa1u_n),grid,xlabel('freq(Hz)'),ylabel('PAa1u (deg)')
subplot(222),plot(freq,db(abs(TFa2u)),freq,db(abs(TFa2u_n))),grid,xlabel('freq(Hz)'),ylabel('TFa2u (dB)'),legend('w/o noise','w/ noise')
subplot(224),plot(freq,PAa2u,freq,PAa2u_n),grid,xlabel('freq(Hz)'),ylabel('PAa2u (deg)')
figure
subplot(221),plot(freq,db(abs(TFa3u)),freq,db(abs(TFa3u_n))),grid,xlabel('freq(Hz)'),ylabel('TFa3u (dB)'),legend('w/o noise','w/ noise')
subplot(223),plot(freq,PAa3u,freq,PAa3u_n),grid,xlabel('freq(Hz)'),ylabel('PAa3u (deg)')

% save the data
freq=freq';
save pexdata TFd1w TFa1w TFa2w TFa3w TFfw TFd1w_n TFa1w_n TFa2w_n TFa3w_n TFfw_n ...
    PAd1w PAa1w PAa2w PAa3w PAfw PAd1w_n PAa1w_n PAa2w_n PAa3w_n PAfw_n ...
    TFd1u TFa1u TFa2u TFa3u TFfu TFd1u_n TFa1u_n TFa2u_n TFa3u_n TFfu_n ...
    PAd1u PAa1u PAa2u PAa3u PAfu PAd1u_n PAa1u_n PAa2u_n PAa3u_n PAfu_n freq
fprintf('Get transfer function data from analytical model using frequency response function ...Done\n');
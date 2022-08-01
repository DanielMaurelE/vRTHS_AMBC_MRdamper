function [H1,coh,fuu] = H1est(y,u,w,nfft,fs)
% Entradas:
%     y    : vector de la señal de salida
%     u    : vector de la señal de entrada
%     w    : tamaño de las ventanas
%     nfft : número de puntos de la FFT
%     fs   : frecuencia de sampleo

% Salidas:
%     H1   : Valor del estimador H1 en número complejo
%     coh  : Valor de la coherencia
%     fuu  : vector de las frecuencias

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[Syu,fyu] = cpsd(u,y,w,[],nfft,fs,'onesided');
[Suu,fuu] = cpsd(u,u,w,[],nfft,fs,'onesided'); 
[Syy,fyy] = cpsd(y,y,w,[],nfft,fs,'onesided');

H1 = Syu./Suu; % dB
H2 = Syy./Syu; % dB
R = corrcoef(H1,H2);
R1 = corrcoef(H1);
R2 = corrcoef(H2);
coh = abs(Syu).^2./(Suu.*Syy);


% Compare results
% coded by K. S. Park at UIUC

clear all; close all; clc;
fprintf('Compare results of analytical model and esitmated model ...\n');
% Compare the transfer function
load ABCD_ESTI_RED
% compare full- and reduced-order model of MIMO obtained from two SIMO (MIMO 1 system)
figure,bode(At,Bt,Ct,Dt,1),hold,bode(Aer1,Ber1,Cer1,Der1,1,'r') % earthquake excitaion
title('MIMO 1 - from excitation input');
legend('Full-order model','Reduced-order model');
figure,bode(At,Bt,Ct,Dt,2),hold,bode(Aer1,Ber1,Cer1,Der1,2,'r') % control input
title('MIMO 1 - from control input');
legend('Full-order model','Reduced-order model');
% compare full- and reduced-order model of MIMO obtained from one MIMO (MIMO 2 system)
figure,bode(At2,Bt2,Ct2,Dt2,1),hold,bode(Aer2,Ber2,Cer2,Der2,1,'r') % earthquake excitation
title('MIMO 2 - from excitation input');
legend('Full-order model','Reduced-order model');
figure,bode(At2,Bt2,Ct2,Dt2,2),hold,bode(Aer2,Ber2,Cer2,Der2,2,'r') % control input
title('MIMO 2 - from control input');
legend('Full-order model','Reduced-order model');
% compare pz plane for MIMO 1 system
[z,p]=ss2zp(At,Bt,Ct,Dt,1);
[zr,pr]=ss2zp(Aer1,Ber1,Cer1,Der1,1);
figure,plot(real(p),imag(p),'rx',real(pr),imag(pr),'bo'),legend('full-order','reduced-order'),title('MIMO 1 from w - Pole');
xlabel('Imag'),ylabel('Real');
figure; 
for i=1:4
    subplot(2,2,i),plot(real(z(:,i)),imag(z(:,i)),'rx',real(zr(:,i)),imag(zr(:,i)),'bo'),legend('full-order','reduced-order');
    xlabel('Imag'),ylabel('Real'), title('MIMO 1 from w - Zeros');
end
[z,p]=ss2zp(At,Bt,Ct,Dt,2);
[zr,pr]=ss2zp(Aer1,Ber1,Cer1,Der1,2);
figure,plot(real(p),imag(p),'rx',real(pr),imag(pr),'bo'),legend('full-order','reduced-order'),title('MIMO 1 from u - Pole')
figure;
for i=1:4
    subplot(2,2,i),plot(real(z(:,i)),imag(z(:,i)),'rx',real(zr(:,i)),imag(zr(:,i)),'bo'),legend('full-order','reduced-order'),title('zero')
    xlabel('Imag'),ylabel('Real'), title('MIMO 1 from u - Zeros');    
end

% compare pz plane for MIMO 2 system
[z,p]=ss2zp(At2,Bt2,Ct2,Dt2,1);
[zr,pr]=ss2zp(Aer2,Ber2,Cer2,Der2,1);
figure,plot(real(p),imag(p),'rx',real(pr),imag(pr),'bo'),legend('full-order','reduced-order'),title('MIMO 2 from w - Pole');
figure;
for i=1:4
    subplot(2,2,i),plot(real(z(:,i)),imag(z(:,i)),'rx',real(zr(:,i)),imag(zr(:,i)),'bo'),legend('full-order','reduced-order');
    xlabel('Imag'),ylabel('Real'), title('MIMO 2 from w - Zeros');
end
[z,p]=ss2zp(At2,Bt2,Ct2,Dt2,2);
[zr,pr]=ss2zp(Aer2,Ber2,Cer2,Der2,2);
figure,plot(real(p),imag(p),'rx',real(pr),imag(pr),'bo'),legend('full-order','reduced-order'),title('MIMO 2 from u - Pole');
figure
for i=1:4
    subplot(2,2,i),plot(real(z(:,i)),imag(z(:,i)),'rx',real(zr(:,i)),imag(zr(:,i)),'bo'),legend('full-order','reduced-order');
    xlabel('Imag'),ylabel('Real'), title('MIMO 2 from u - Zeros');
end

% compare MIMO 1 and MIMO 2 system (full-order model)
[z1,p1]=ss2zp(At,Bt,Ct,Dt,1);
[z2,p2]=ss2zp(At2,Bt2,Ct2,Dt2,1);
figure,plot(real(p1),imag(p1),'rx',real(p2),imag(p2),'bo'),legend('MIMO 1','MIMO 2'),title('full-order model from w - Pole')
figure
for i=1:4
    subplot(2,2,i),plot(real(z1(:,i)),imag(z1(:,i)),'rx',real(z2(:,i)),imag(z2(:,i)),'bo'),legend('MIMO 1','MIMO 2'),title('full-order model from w - Zeros');
end
[z1,p1]=ss2zp(At,Bt,Ct,Dt,2);
[z2,p2]=ss2zp(At2,Bt2,Ct2,Dt2,2);
figure,plot(real(p1),imag(p1),'rx',real(p2),imag(p2),'bo'),legend('MIMO 1','MIMO 2'),title('full-order model from u- Pole')
figure
for i=1:4
    subplot(2,2,i),plot(real(z1(:,i)),imag(z1(:,i)),'rx',real(z2(:,i)),imag(z2(:,i)),'bo'),legend('MIMO 1','MIMO 2'),title('full-order model from u- Zeros')
end

% compare eigenvalue
load ABCD
eig_anal = eig(A1);
eig_mimo1 = eig(Aer1);
eig_mimo2 = eig(Aer2);
tmp = [eig_anal zeros(7,1) eig_mimo1 zeros(7,1) eig_mimo2]
fprintf('==== Analtytical Model ====        ==== Estimated Model 1 ====           ==== Esitmated Model 2 ====\n');

fprintf('\nCompare results of analytical model and esitmated model ...Done\n');

% Reduce the model
% coded by K. S. Park at UIUC

clear all; close all; clc;
fprintf('Reduce model from full-order model to remove repeated poles ...\n');
load ABCD_ESTI              % load estimated model
% Combine two SIMO system as one MIMO system
At=daug(Aew,Aeu);
Bt=daug(Bew,Beu);
Ct=[Cew Ceu];
Dt=[Dew Deu];

% data from MIMO system
At2=Ae;
Bt2=Be;
Ct2=Ce;
Dt2=De;

% converting the model to canonical modal form
[At,Bt,Ct,Dt,t]=canon(At,Bt,Ct,Dt,'modal');
% balancing the model
[Ab,Bb,Cb,g,t]=balreal(At,Bt,Ct);
plot(g),xlabel('state'),ylabel('Hankel singular value'),grid
N=input('Number of State to be remained ? ');
elim=[N+1:size(At,1)];
%elim=(g<1e-1);
% form reduced model
[Aer1,Ber1,Cer1,Der1]=modred(Ab,Bb,Cb,Dt,elim);

% converting the model to canonical modal form
[At2,Bt2,Ct2,Dt2,t]=canon(At2,Bt2,Ct2,Dt2,'modal');
% balancing the model
[Ab,Bb,Cb,g,t]=balreal(At2,Bt2,Ct2);
plot(g),xlabel('state'),ylabel('Hankel singular value'),grid
N=input('Number of State to be remained ? ');
elim=[N+1:size(At2,1)];
%elim=(g<1e-1);
% form reduced model
[Aer2,Ber2,Cer2,Der2]=modred(Ab,Bb,Cb,Dt2,elim);

save ABCD_ESTI_RED At Bt Ct Dt At2 Bt2 Ct2 Dt2 Aer1 Ber1 Cer1 Der1 Aer2 Ber2 Cer2 Der2
fprintf('Reduce model from full-order model to remove repeated poles ...Done\n');
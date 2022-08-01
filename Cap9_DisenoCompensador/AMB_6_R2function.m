% Function for evaluation of multiple simulations with one matrix of
% adaptive gains
function [R2,dR2,maxJ2] = AMB_6_R2function(x,escalas,frefcals,drefcals,realizations,b1cals,b2cals,b3cals)
%Inputs:
% x       :  1x4 vector x=log10(diag(Gamma))
%escalas  :  1x2 vector earthquake scale [minScale maxScale]
%frefcalas :  1x2 vector SDOF numerical structure natural frequency range
%[fmin fmax] (Hz)
%drefclas  :  1x2 vector SDOF damping ratio range [dmin dmax]
%mecals  :  1x2 vector experimental mass range [m_min m_max]
%cecals  :  1x2 vector experimental viscous damping range [cmin cmax]
%kecals  :  1x2 vector experimental stiffness range [kmin kmax]
%realizations  : Number of simulations


%Outputs:
%R2 : mean(J2)
%dR2 : std(J2)
%maxJ2  : max(J2)

%Function:
adaptivegain=diag(10.^x);   %define adaptive gains matrix
assignin('base','adaptivegain',adaptivegain);
    
    %Define scales, SDOF properties and exp. sub. properties for the N
    %realizations
    latin = lhsdesign(realizations,6); %latin hypercube sampling
    escalapar = escalas(1)+(escalas(2)-escalas(1))*latin(:,1);
    drefcalpar = drefcals(1)+(drefcals(2)-drefcals(1))*latin(:,2);
    wnrefcalpar = 2*pi*frefcals(1)+2*pi*(frefcals(2)-frefcals(1))*latin(:,3);
    b1par = b1cals(1)+(b1cals(2)-b1cals(1))*latin(:,4);
    b2par = b2cals(1)+(b2cals(2)-b2cals(1))*latin(:,5);
    b3par = b3cals(1)+(b3cals(2)-b3cals(1))*latin(:,6);

% Parallel simulations
J2s = zeros(realizations,1);
for i=1:realizations
    in(i) = Simulink.SimulationInput('AMB_9_calibration_simulation');
    in(i) = in(i).setVariable('escala',escalapar(i));
    in(i) = in(i).setVariable('drefcal',drefcalpar(i));
    in(i) = in(i).setVariable('wnrefcal',wnrefcalpar(i));
    
    den = [b3par(i),b2par(i),b1par(i)];
    G = tf(1,den);
    [numcalpar,dencalpar]=tfdata(G,'v');
    
    
    in(i) = in(i).setVariable('numcal',numcalpar);
    in(i) = in(i).setVariable('dencal',dencalpar);
end

    out=parsim(in,'ShowSimulationManager','off','ShowProgress','off','TransferBaseWorkspaceVariables','on');

%J2 evaluation
for i = 1:realizations
        simOut = out(i);
        x_t = simOut.get('x_t').Data;
        x_m = simOut.get('x_m').Data;
        J2s(i)=rms(x_t-x_m)/rms(x_t)*100;
end
%Results
R2 = mean(J2s);
dR2 = std(J2s);
maxJ2 = max(J2s);
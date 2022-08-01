function aux = sincronize(Datos,Ts)
Data(:,3) = Datos(:,2)*9806.65; % Fuerza medida en N
Data(:,2) = Datos(:,1)*10; % desp medido en mm
Data(:,1) = Datos(:,3)*10; % desp comandado en mm

%Noise filter 
n = 4;   %order
fc = 200; %cutoff frequency
[b,a] = butter(n,fc*(Ts/2)); %discrete filter coefficients
Data = filter(b,a,Data);

Data = detrend(Data,0);
Data = Data-mean(Data(end-40:end,:));
pos = find(abs(Data(:,1))>0.4);

D = Data(pos(1)-1000:pos(end)+1000,:);


x_c = D(1:end,1);
x_m = D(1:end,2);
F = D(1:end,3);
aux = [x_c,x_m,F];
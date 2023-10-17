clearvars;
clc;
L = 64; %number of digital samples per data bit
Fs = 10*L; %Sampling frequency
voltageLevel = 5;
data = rand(1000,1)>0.5;
clk = mod(0:2*numel(data)-1,2).';
ami = 1*data;
previousBit = 0; %AMI encoding
for idx=1:numel(data)
if (ami(idx)==1) && (previousBit==0)
ami(idx)= voltageLevel;
previousBit=1;
end
if (ami(idx)==1) && (previousBit==1)
ami(idx)= -voltageLevel;
previousBit = 0;
end
end
clk_sequence=reshape(repmat(clk,1,L).',1,length(clk)*L);
data_sequence=reshape(repmat(data,1,2*L).',1,length(data)*2*L);
unipolar_nrz_l = voltageLevel*data_sequence;
nrz_encoded = voltageLevel*(2*data_sequence - 1);
unipolar_rz = voltageLevel*and(data_sequence,not(clk_sequence));
ami_sequence = reshape(repmat(ami,1,2*L).',1,length(ami)*2*L);
manchester_encoded = voltageLevel*(2*xor(data_sequence,clk_sequence)-1);
figure(1);
subplot(4,1,1);
plot(clk_sequence(1:800),"LineWidth",3);
title('Clock');
grid on
subplot(4,1,2);
plot(data_sequence(1:800),"LineWidth",3);
title('Data');
grid on
subplot(4,1,3);
plot(unipolar_nrz_l(1:800),"LineWidth",3);
title('Unipolar non-return-to-zero level');
grid on;
subplot(4,1,4);
plot(nrz_encoded(1:800),"LineWidth",3);
title('Bipolar Non-return-to-zero level');
grid on;
figure(2);
subplot(4,1,1);
plot(clk_sequence(1:800),"LineWidth",3);
title('Clock');
grid on;
subplot(4,1,2);
plot(data_sequence(1:800),"LineWidth",3);
title('Data');
grid on;
subplot(4,1,3);
plot(ami_sequence(1:800),"LineWidth",3);
title('Alternate Mark Inversion (AMI)');
grid on;
subplot(4,1,4);
plot(manchester_encoded(1:800),"LineWidth",3);
title('Manchester Scheme');
grid on;
Rb=1;
Tb=1/Rb;
f=0:0.025*Rb:2*Rb;
x=f*Tb;
P1 = (Tb)*(sinc(x).*sinc(x));
figure(1);
plot(f,P1,'r');
grid on;
box on;
xlabel('frequency as a multiple of Bitrate(fRb)---->');
ylabel('Power Spectral Density ---->');
title('PSD for Polar Signal');
hold on;
P2 = (Tb)*(sinc(x/2)).*(sinc(x/2)).*(sin(pi*x)).*(sin(pi*x));
plot(f,P2,'m');
grid on;
box on;
xlabel('frequency as a multiple of Bitrate(fRb)---->');
ylabel('Power Spectral Density ---->');
title('PSD for Polar Signal');
hold on;
P3 = (Tb/4)*(sinc(x).*sinc(x)) + (1/4)*dirac(f);
plot(f,P3,'g');
grid on;
box on;
xlabel('frequency as a multiple of Bitrate(fRb)---->');
ylabel('Power Spectral Density ---->');
title('PSD for Polar Signal');
hold on;
P4 = Tb*(sinc(x/2)).*(sinc(x/2)).*(sin(pi*x/2)).*(sin(pi*x/2));
plot(f,P4,'b');
grid on;
box on;
xlabel('frequency as a multiple of Bitrate(fRb)---->');
ylabel('Normalised Power Spectral Density ---->');
title('PSD for Different Line Coding Schemes');
hold on;
legend('Polar','Bipolar','Unipolar','Manchester');
clear;
clc; %Signal Generation
b = input('Enter the Bit stream \n ');
n = length(b);
t = 0:.01:n;
x = 1:1:(n+1)*100;
f = 20;
sint = sin(2*pi*f*t);
bw1 = zeros(1, length(t));
bw2 = zeros(1, length(t));
for i = 1:n
for j = i:.1:i+1
bw1(x(i*100:(i+1)*100)) = b(i);
bw2(x(i*100:(i+1)*100)) = b(i);
end
f = f + 12;
sint = sin(2*pi*f*t);
end
bw1 = bw1(100:end);
bw2 = bw2(100:end);
figure;
subplot(4,1,1)
plot(t,bw1)
grid on;
axis([0 n -2 +2])
title('Modulating Signal')
subplot(4,1,2)
st1 = bw1.*sint;
plot(t, st1)
grid on;
axis([0 n -2 +2])
title('Modulated Signal ')
for i = 1:n
if (b(i) == 0)
b_p(i) = -1;
else
b_p(i) = 1;
end
for j = i:.1:i+1
bw2(x(i*100:(i+1)*100)) = b_p(i);
end
f = f + 12;
sint = sin(2*pi*f*t);
end
bw2 = bw2(100:end);
wo = 2*(2*pi*t);
W = 1*(2*pi*t);
sinHt = sin(wo+W);
sinLt = sin(wo-W);
st2 = sin(wo+(bw2).*W);
subplot(4,1,3)
plot(t, st2)
grid on;
axis([0 n -2 +2])
title('Modulated Signal BFSK')
for i = 1:n
if (b(i) == 0)
b_p2(i) = -1;
else
b_p2(i) = 1;
end
for j = i:.1:i+1
bw2(x(i*100:(i+1)*100)) = b_p2(i);
end
f = f + 12;
sint = sin(2*pi*f*t);
end
bw2 = bw2(100:end);
st3 = bw2.*sint;
subplot(4,1,4)
plot(t, st3)
grid on;
axis([0 n -2 +2])
title('Modulated Signal BPSK')
Eb_N0dB=0:20;
EbN0=10.^(Eb_N0dB/10);
pe_bfsk=0.5*erfc(sqrt(EbN0/2));
figure;
subplot(3,1,1);
semilogy(Eb_N0dB,pe_bfsk)
xlabel('Eb/N0(dB)');
ylabel('BER');
title('BFSK BER RATE GRAPH - THEORETICAL ')
axis([0 20 0.000001 1])
grid on;
Eb_N0dB=0:20;
EbN0=10.^(Eb_N0dB/10);
pe_bpsk=0.5*erfc(sqrt(EbN0));
subplot(3,1,2);
semilogy(Eb_N0dB,pe_bpsk)
title('BPSK BER RATE GRAPH - THEORITICAL')
xlabel('Eb/N0(dB)');
ylabel('BER');
axis([0 20 0.000001 1])
grid on;
Eb_N0dB=0:20;
EbN0=10.^(Eb_N0dB/10);
pe_bask=0.5*erfc(sqrt(EbN0/4));
subplot(3,1,3);
semilogy(Eb_N0dB,pe_bask)
title('BASK BER RATE GRAPH - THEORITICAL')
xlabel('Eb/N0(dB)');
ylabel('BER');
axis([0 20 0.000001 1])
grid on;
num_bit = 10000;
BER_iter = 20;
Eb=1;
SNRdB=0:0.1:10;
SNR=10.^(SNRdB/10);
for count=1:length(SNR)
avgError=0;
No=Eb/SNR(count);
for run_time=1:BER_iter
Error=0;
data = randi([0 1],1,num_bit);
Y = awgn(complex(data),SNRdB(count));
for k=1:num_bit
if ((Y(k)>0.5 && data(k)==0)||(Y(k)<0.5 && data(k)==1))
Error=Error+1;
end
end
Error=Error/num_bit;
avgError=avgError+Error;
end
BER_sim(count)=avgError/BER_iter;
end
figure (4)
semilogy(SNRdB,BER_sim,'g','linewidth',2.5);
grid on;
hold on;
BER_th = (1/2)*erfc(0.5*sqrt(SNR));
semilogy(SNRdB,BER_th,'r','linewidth',2.5);
grid on;
hold on;
title(' Curve for Bit Error Rate verses SNR for ASK modulation');
xlabel(' SNR(dB)');
ylabel('BER');
legend('Simulation','Theoretical');
num_bit = 1000;
BER_iter = 20;
Eb=1;
SNRdB=0:0.2:10;
SNR=10.^(SNRdB/10);
for count=1:length(SNR)
avgError=0;
No=Eb/SNR(count);
for run_time=1:BER_iter
Error=0;
data = randi([0 1],1,num_bit);
s=data+1i*(~data);
Nimg = sqrt(No/2)*randn(1,num_bit);
Nreal = sqrt(No/2)*randn(1,num_bit);
N = Nimg+1i*Nreal;
Y = s+N;
for k=1:num_bit
Z(k)=real(Y(k))-imag(Y(k));
if ((Z(k)>0 && data(k)==0)||(Z(k)<0 && data(k)==1))
Error=Error+1;
end
end
Error=Error/num_bit;
avgError=avgError+Error;
end
BER_sim(count)=avgError/BER_iter;
end
figure (5)
semilogy(SNRdB,BER_sim,'g','linewidth',2.5);
grid on;
hold on;
BER_th=(1/2)*erfc(sqrt(SNR/2));
semilogy(SNRdB,BER_th,'r','linewidth',2.5);
grid on;
hold on;
title(' Curve for Bit Error Rate verses SNR for BFSK modulation');
xlabel(' SNR(dB)');
ylabel('BER');
legend('Simulation','Theoretical') ;
%BPSK MODULATION
num_bit = 1000;
BER_iter = 20;
Eb=1;
SNRdB=0:0.2:10;
SNR=10.^(SNRdB/10);
for count=1:length(SNR)
avgError=0;
No=Eb/SNR(count);
for run_time=1:BER_iter
Error=0;
data = randi([0 1],1,num_bit);
s=2*data-1;
N = sqrt(No/2)*randn(1,num_bit);
Y = s+N;
for k=1:num_bit
if ((Y(k)>0 && data(k)==0)||(Y(k)<0 && data(k)==1))
Error=Error+1;
end
end
Error=Error/num_bit;
avgError=avgError+Error;
end
BER_sim(count)=avgError/BER_iter;
end
figure (6)
semilogy(SNRdB,BER_sim,'g','linewidth',2.5);
grid on;
hold on;
BER_th=(1/2)*erfc(sqrt(SNR));
semilogy(SNRdB,BER_th,'r','linewidth',2.5);
grid on;
hold on;
title(' Curve for Bit Error Rate verses SNR for BPSK modulation');
xlabel(' SNR(dB)');
ylabel('BER');
legend('Simulation','Theoretical');
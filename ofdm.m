% OFDM transmitter

% Define parameters
numCarriers = 64;
cpLength = 16;
numSymbols = 10;

% Generate random data stream
data = randi([0, 1], numSymbols, numCarriers);

% Convert data to time domain
timeDomainSignal = ifft(data, numCarriers, 2);

% Add cyclic prefix
timeDomainSignalWithCP = [timeDomainSignal(:, end - cpLength + 1:end), timeDomainSignal];

% Plot time domain signal with CP
figure;
plot(0:length(timeDomainSignalWithCP)-1,real(timeDomainSignalWithCP(1,:)));
title('Time Domain OFDM Symbol');
xlabel('Sample Index');
ylabel('Amplitude');

% Generate frequency domain signal
freqDomainSignal = fft(timeDomainSignalWithCP, numCarriers, 2);

% Plot frequency domain signal
figure;
plot(0:numCarriers-1, abs(freqDomainSignal(1, :)));
title('Frequency Domain OFDM Symbol');
xlabel('Subcarrier Index');
ylabel('Amplitude');

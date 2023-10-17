% Rake Receiver Simulation (First 20 values)
% Parameters
numBits = 20; % Number of transmitted bits (adjusted for plotting)
snrdB = 10; % Signal-to-noise ratio in dB
numPaths = 3; % Number of multipath components
delayVector = [0, 1, 2]; % Delays of the multipath components
numTaps = max(delayVector) + 1; % Number of taps in the rake receiver
% Generate random bits for transmission
transmittedBits = randi([0, 1], 1, numBits);
% BPSK Modulation
modulatedSymbols = 2 * transmittedBits - 1;
% Add Gaussian noise to the signal
noisePower = 10^(-snrdB / 10);
receivedSignal = awgn(modulatedSymbols, snrdB, 'measured');
% Rake Receiver
rakeOutput = zeros(1, numBits);
for i = 1:numBits
% Create delayed versions of the received signal
delayedSignals = zeros(1, numPaths);
for j = 1:numPaths
if i - delayVector(j) > 0
delayedSignals(j) = receivedSignal(i - delayVector(j));
end
end
% Rake combining - take the sign of the sum
rakeOutput(i) = sign(sum(delayedSignals));
end
disp(rakeOutput)
% Decision making (hard decision)
receivedBits = (rakeOutput > 0);
% Display results
disp(['Transmitted Bits: ', num2str(transmittedBits)]);
disp(['Received Bits: ', num2str(receivedBits)]);
% Plot transmitted and received bits for the first 20 values
figure(1)
subplot(2,1,1)
stem(1:numBits, transmittedBits, 'b', 'LineWidth', 1.5);
title('Transmitted Bits (First 20 values)');
xlabel('Bit Index');
ylabel('Bit Value');
subplot(2,1,2)
stem(1:numBits, receivedBits, 'r', 'LineWidth', 1.5);
title('Transmitted Bits (First 20 values)');
xlabel('Bit Index');
ylabel('Bit Value');
hold off;
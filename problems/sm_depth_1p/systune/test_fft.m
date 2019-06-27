Fs = 10;            % Sampling frequency                    
T = 1/Fs;             % Sampling period       
L = 150000;             % Length of signal
t = (0:L-1)*T;        % Time vector
K = 50;
X = K*randn(size(t));
figure;plot(10*t(1:50),X(1:50));
title('Signal Corrupted with Zero-Mean Random Noise');
xlabel('t (milliseconds)');
ylabel('X(t)');

Y = fft(X);
P2 = abs(Y/L);
P1 = P2(1:L/2+1);
P1(2:end-1) = 2*P1(2:end-1);
f = Fs*(0:(L/2))/L;
figure;plot(f,P1);
title('Single-Sided Amplitude Spectrum of X(t)');
xlabel('f (Hz)');
ylabel('|P1(f)|');
mean(P1)
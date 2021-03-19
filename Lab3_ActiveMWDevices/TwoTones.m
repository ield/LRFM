%% Dos tonos
% SCRIPT to create a configuration file for MXG generator that generates
% two tones for the calculation of third order intermodulation products.

clear; clc; close all;

fs = 125*10^6;    % MXG sampling frequency. Max is 125 MSPS.
f_tone=5*10^6;    % Separation of each tone from the carrier

NumSamples=625;   % Number of samples of the sampled modulating signal 
n=0:NumSamples-1; % Sample indices
t=n/fs;           % Time vector

signal=exp(1j*2*pi*f_tone*t+1j*2*pi*rand())+exp(1j*2*pi*(-f_tone)*t+1j*2*pi*rand());
signalFFT=fft(signal.*hann(NumSamples)',1024);          % Hanning windowing

figure('Color',[1 1 1]);

subplot(2,1,1)
plot(cat(2,((0:511)*fs/1024*10^-6)-fs/2*10^-6,(0:511)*fs/1024*10^-6),...
    cat(2,20*log10(abs(signalFFT(512:end))),20*log10(abs(signalFFT(1:511))))-49.87);
grid; axis([-30 30 -150 0]);
xlabel('f(MHz)')
ylabel('dBm');

subplot(2,1,2)
plot(t,real(signal)); grid
hold on;
plot(t,imag(signal));
legend('I-Component','Q-Component')

vector_a_agilent(signal,'2TONES'); % Function to generate the MXG readable file

set(gcf,'position',[100,100,400,300]);
set(gca,'FontSize',11, 'fontname','Times New Roman');

save_path = '../../Lab3_ActiveMWDevices/Images/';
saveas(gca, [save_path, '2_tones'],'epsc');
%% MULTITONO
% SCRIPT to create a file for MXG generator that generates N equally spaced
% tones for the calculation of M-IMR and ACPR

clear; clc; close all

fs=125*10^6;        % MXG sampling frequency. Max is 125 MSPS.
delta_f=1*10^6;     % Tone-to-tone separation
NumPortadoras=29;   % Number of tones (odd)

NumMuestras=625;    % Number of samples of the sampled modulating signal 
n=0:NumMuestras-1;  % Sample indices
t=n/fs;             % Time vector

signal=ones(1,625); % The central tone is at DC before modulation

numTones=floor((NumPortadoras-1)/2);

tone_to_delete = 3;
% tone_to_delete = 20;


for n=1:numTones
    if n == tone_to_delete
        signal=signal+exp(1j*2*pi*n*delta_f*t+1j*2*pi*rand());
    else
        signal=signal+exp(1j*2*pi*n*delta_f*t+1j*2*pi*rand())+...
            exp(-1j*2*pi*n*delta_f*t+1j*2*pi*rand());
    end
end
numTonesWithCarrier=2*numTones + 1;

signalFFT=fft(signal.*hann(NumMuestras)',1024*10);

figure('Color',[1 1 1]);

subplot(2,1,1)
plot(cat(2,((0:5119)*fs/10240*10^-6)-fs/2*10^-6,(0:5119)*fs/10240*10^-6),...
    cat(2,20*log10(abs(signalFFT(5120:end))),20*log10(abs(signalFFT(1:5119))))-49.87);
axis([-20 20 -140 5])
grid;
xlabel('f(MHz)')
ylabel('dBm');

vector_a_agilent(signal,['MULTITONE1_',num2str(numTonesWithCarrier),'TONES',num2str(delta_f/1e6),'MHZ'])

subplot(2,1,2)
plot(t,real(signal)); hold on
plot(t,imag(signal));
legend('I-Component','Q-Component')
axis tight
grid

set(gcf,'position',[100,100,400,300]);
set(gca,'FontSize',11, 'fontname','Times New Roman');

save_path = '../../Lab3_ActiveMWDevices/Images/';
% saveas(gca, [save_path, 'multitone_notch'],'epsc');
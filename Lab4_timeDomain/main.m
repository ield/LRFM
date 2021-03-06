close all; clc; clear;
save_path = '../../Lab4_timeDomain/Images/';
% Obtain the data from the measurement file
[s11, f_original] = getSParametersFromFile_Reflection('filter_wrong.dat');

% Obtain the value at f = 0 by extrapollating
% f_complete contains f = 0
f_complete = [0; f_original];
s11_complete = interp1(f_original, s11, f_complete, 'linear', 'extrap');

% Initialize the fft time domain and frequency domain axis
[t, f] = initializeFT(length(f_complete)-1,f_complete(2));

figure('Color',[1 1 1]);
set(gcf,'position',[100,100,400,300]);

plot(f*1e9, 20*log10(abs(s11_complete)));
xlim([f(1) f(end)]*1e9);
xlabel('Frequency (GHz)');
ylabel('s_{11} (dB)');
set(gca,'FontSize',11, 'fontname','Times New Roman');
saveas(gca, [save_path, 'measured_freq'],'epsc');

% Apply windowing to the signal to reduce the effects on
% truncation/alliasing. It is started with a kaiser window of order k
k = 6;
window_kaiser = kaiserwindowuni(length(f)-1, k);
s11_windowed = s11_complete .* window_kaiser;

figure('Color',[1 1 1]);
set(gcf,'position',[100,100,400,300]);

plot(f*1e9, 20*log10(abs(s11_windowed))); hold on;
plot(f*1e9, 20*log10(abs(window_kaiser))); hold on;
plot(f*1e9, 20*log10(abs(s11_complete)), ':');
xlim([f(1) f(end)]*1e9);
xlabel('Frequency (GHz)');
ylabel('s_{11} (dB)');
legend('Windowed s_{11}', 'Kaiser Window', 'Measured signal', ...
    'location', 'southwest');
set(gca,'FontSize',11, 'fontname','Times New Roman');
saveas(gca, [save_path, 'windowed'],'epsc');


% Obtain the signal in the time domain
s11_time = rifftuni(s11_windowed);

figure('Color',[1 1 1]);
set(gcf,'position',[100,100,800,300]); 

subplot(1, 2, 1);
plot(t*1e9, s11_time);
xlim([t(1) t(end)]*1e9);
xlabel('Time (ns)');
set(gca,'FontSize',11, 'fontname','Times New Roman');

subplot(1, 2, 2);
plot(t*1e9, s11_time)
xlabel('Time (ns)');
xlim([0 20]);
set(gca,'FontSize',11, 'fontname','Times New Roman');

saveas(gca, [save_path, 'measured_time'],'epsc');

% It is seen that the perturbations of the connector end in t = 1ns.
% Therefore, the gate should elliminate these frequencies.
% Find out the time in which the filter should be
length_cable = 42e-2;   %(m)
c_cable = 0.6*physconst('LightSpeed');
t_filter = length_cable/c_cable;
fprintf('t = %f ns', t_filter*1e9);

% It is determined the gate used: the effects of the connector must be
% corrected before and after the filter
t_ini = 5e-9;
t_fin = 10.5e-9;
tram = 0.5e-9;
gate_filter = gate(t, t_ini, t_fin, tram);
s11_time_gated = s11_time .* gate_filter;

figure('Color',[1 1 1]);
set(gcf,'position',[100,100,800,300]); 

subplot(1, 2, 1);
plot(t*1e9, s11_time_gated/max(s11_time)); hold on;
plot(t*1e9, gate_filter); hold on;
plot(t*1e9, s11_time/max(s11_time), ':k'); hold on;
xlim([t(1) t(end)]*1e9);
legend('Normalized gated, windowed signal', 'Gate', ...
    'Normalized windowed signal', 'location', 'southwest');

xlabel('Time (ns)');
set(gca,'FontSize',11, 'fontname','Times New Roman');

subplot(1, 2, 2);
plot(t*1e9, s11_time_gated/max(s11_time)); hold on;
plot(t*1e9, gate_filter); hold on;
plot(t*1e9, s11_time/max(s11_time), ':k'); hold on;

xlim([0 25]);
xlabel('Time (ns)');
set(gca,'FontSize',11, 'fontname','Times New Roman');

saveas(gca, [save_path, 'gated'],'epsc');

% Obtain the corrected signal in the frequency domain
s11_corrected_windowed = fftuni(s11_time_gated);

figure('Color',[1 1 1]);
set(gcf,'position',[100,100,400,300]);

plot(f*1e9, 20*log10(abs(s11_corrected_windowed))); hold on;
plot(f*1e9, 20*log10(abs(s11_complete)), ':');
xlim([f(1) f(end)]*1e9);
xlabel('Frequency (GHz)');
ylabel('s_{11} (dB)');
legend('Windowed, corrected s_{11}', 'Measured signal', ...
    'location', 'southwest');
set(gca,'FontSize',11, 'fontname','Times New Roman');
saveas(gca, [save_path, 'windowed_corrected'],'epsc');


% Undo the windowing performed
s11_corrected = s11_corrected_windowed./window_kaiser;

figure('Color',[1 1 1]);
set(gcf,'position',[100,100,400,300]);

plot(f*1e9, 20*log10(abs(s11_corrected))); hold on;
plot(f*1e9, 20*log10(abs(s11_complete)), ':');
xlim([f(1) f(end)]*1e9);
xlabel('Frequency (GHz)');
ylabel('s_{11} (dB)');
legend('Corrected s_{11}', 'Measured signal', ...
    'location', 'southwest');
set(gca,'FontSize',11, 'fontname','Times New Roman');
saveas(gca, [save_path, 'corrected'],'epsc');



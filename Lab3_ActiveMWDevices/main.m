%% Author: ield
% Subject: LRFM MSTC ETSIT UPM
%% Description of the code
% This code gathers the work done in order to characterize an active device
% 1. S-parameters of the characterization element
% 2. S-parameters of the amplifier.
clear;
clc;
close all;
save_path = '../../Lab3_ActiveMWDevices/Images/';

%% 1. S-parameters of the characterization element
% Obtain the S-parameters
[f, s11_mod, s11_phase, s21_mod, s21_phase, s12_mod, s12_phase, ...
    s22_mod, s22_phase] = getSParametersFromFile_s2p('calibration_att.s2p');

% Plot the S-parameters
plotS_Parameters_2port(f, s11_mod, s11_phase, s21_mod, s21_phase, ...
    s12_mod, s12_phase, s22_mod, s22_phase, -60, -5, ...
    save_path, 'att_calibration');

%% 2. S-parameters of the amplifier.
% Obtain the S-parameters
[f, s11_mod, s11_phase, s21_mod, s21_phase, s12_mod, s12_phase, ...
    s22_mod, s22_phase] = getSParametersFromFile_s2p('s_par_amplifier.s2p');

% Plot the S-parameters
plotS_Parameters_2port(f, s11_mod, s11_phase, s21_mod, s21_phase, ...
    s12_mod, s12_phase, s22_mod, s22_phase, -30, 20, ...
    save_path, 's_par_amplifier');

%% 3. P_1dB of the amplifier.
% Obtain the S-parameters
[p, ~, ~, s21_mod, s21_phase, ~, ~, ~, ~] = ...
    getSParametersFromFile_s2p('s21_power_db.s2p');

% Plot the S-parameters
plot_power_sweep(p, s21_mod, s21_phase, save_path, 'p1db');

%% 4. Input power drift of the VNA
% Obtain the measured values
[measured_results] = readtable('OpcionalBarridoPotencia.xlsx');
% Convert from table to matrix
measured_results = measured_results{:,:};
% There are results starting at -23dBm. We are only interested in the
% results starting in -20 dBm
P_vna = measured_results(4:end,1);
P_measured = measured_results(4:end,2);

% Plot the drift
% Plot the S-parameters of the device (phase)
figure('Color',[1 1 1]);

plot(P_vna, P_measured); hold on;
xlabel('Ideal power (dBm)');
xlim([P_vna(1) P_vna(end)]);
ylabel('Measured power (dBm)');

set(gcf,'position',[100,100,400,300]);
set(gca,'FontSize',11, 'fontname','Times New Roman');

saveas(gca, [save_path, 'power_drift'],'epsc');

% Plot the S-parameters after the drift. In order to plot them, it is
% necessary to make an interpolation of the measured values so the sizes of
% the arrays are equal
P_measured = interp1(P_vna, P_measured, p);
plot_power_sweep(P_measured, s21_mod, s21_phase, save_path, ...
    'p1db_corrected');

%% 5. Noise factor of the amplifier
% Obtain the measurement values and convert to array
measured_results = readtable('TRACE032.CSV', 'HeaderLines',15); 
measured_results = measured_results{:,:};

f = measured_results(:,1);               % Frequency axis
noise_floor_sa = measured_results(:,2);  % Amplifier off,    source off
noise_amplifier = measured_results(:,3); % Amplifier on,     source off
noise_source = measured_results(:,4);    % Amplifier on,     source on

% Plot the three noise sources to see that the noise of the analyzer
% is neglectable
figure('Color',[1 1 1]);

plot(f/1e9, noise_floor_sa); hold on;
plot(f/1e9, noise_amplifier); hold on;
plot(f/1e9, noise_source); hold on;
xlabel('Frequency (GHz)');
xlim([f(1) f(end)]/1e9);
ylabel('Power (dBm)');
legend('Amp: off, source: off', 'Amp: on, source: off', ...
    'Amp: on, source: on', 'location', 'best');

set(gcf,'position',[100,100,400,300]);
set(gca,'FontSize',11, 'fontname','Times New Roman');

saveas(gca, [save_path, 'noise_level'],'epsc');

% Calculate the Y factor (slide 7 of Lab 3). For this, it subtracted the 
% noise floor of the spectrum analyzer (even though this is not necessary)
noise_floor_sa = 10.^(noise_floor_sa/10);
noise_amplifier = 10.^(noise_amplifier/10) - noise_floor_sa;
noise_source = 10.^(noise_source/10) - noise_floor_sa;
Y = noise_source./noise_amplifier;

% In order to find F (slide 7 of Lab 3) it is necessary the ENR, which is
% considered only in the frequencies measured
freq_enr = [1 2 3]*1e9;          % Frequencies GHz
ENR = [15.03 14.96 14.82];       % ENR in dB
% The different values of the ENR are interpoled in order to have arrays of
% the same length of the measurements taken
ENR_interp = interp1(freq_enr, ENR, f);
% Then, the ENR is converted from log to natural units
ENR_interp = 10.^(ENR_interp/10);
% F is calculated (slide 7 of Lab 3) and converted to dB
F_cascade = 10*log10(abs(ENR_interp./(Y-1)));

% To extract the noise figure of the amplifier from the noise figure of the
% cascade elements it is necessary the gain of the amplifier
[f_amp, ~, ~, g_amp, ~, ~, ~, ~, ~] = ...
    getSParametersFromFile_s2p('s_par_amplifier.s2p');
% The measurements are not taken in the same axis. Therefore, first, the
% frequency range is reduced to the frequency range of the noise
% measurements
loc_f_amp_min = find(f_amp >= f(1), 1);
loc_f_amp_max = find(f_amp >= f(end), 1);
f_amp = f_amp(loc_f_amp_min:loc_f_amp_max);
g_amp = g_amp(loc_f_amp_min:loc_f_amp_max);
% Once the frequency range is the same, it is necessary to have the same
% number of points. This is done by interpolating
g_amp = interp1(f_amp, g_amp, f);

% In the connectors the noise figure matches the losses
F1 = 0.2;               % Noise figure of connector 1 in dB
f1 = 10.^(F1/10);       % Noise figure of connector 1 in linear
G1 = -0.2;              % Gain of connector 1 in dB
g1 = 10.^(G1/10);       % Gain of connector 1 in linear

G2 = g_amp;             % Gain of the amplifier in dB
g2 = 10..^(G2/10);      % Gain of the amplifier in linear

F3 = 0.2;               % Noise figure of connector 2 in dB
f3 = 10.^(F3/10);       % Noise figure of connector 2 in linear
G3 = -0.2;              % Gain of connector 2 in dB
g3 = 10.^(G3/10);       % Gain of connector 2 in linear

F4 = 5;                 % Noise figure of sa amplifier in dB
f4 = 10.^(F4/10);       % Noise figure of sa amplifier in linear

f_cascade = 10.^(F_cascade/10);         % Noise figure of cascade elements in linear

% Now it is calculated the noise figure of the amplifier alone
f2 = g1*(f_cascade - f1 - (f3-1)./(g1*g2) - (f4-1)./(g1*g2/g3)) + 1;
F2 = 10*log10(abs(f2));      % Noise figure of the amplifier in dB


% It is plotted the noise figure of the cascade elements
figure('Color',[1 1 1]);

plot(f/1e9, F_cascade, ':'); hold on;
plot(f/1e9, F2, ':'); hold on;
plot(f/1e9, lp_filter_data(F_cascade), 'LineWidth', 2); hold on;
plot(f/1e9, lp_filter_data(F2), 'LineWidth', 2); hold on;
xlabel('Frequency (GHz)');
xlim([f(1) f(end)]/1e9);
ylabel('Noise figure (dB)');
legend('All elements', 'Amplifier', 'All elements filtered', ...
    'Amplifier filtered','location', 'northwest');

set(gcf,'position',[100,100,400,300]);
set(gca,'FontSize',11, 'fontname','Times New Roman');

saveas(gca, [save_path, 'noise_figure'],'epsc');

ylim([2 5]);
saveas(gca, [save_path, 'noise_figure_zoom'],'epsc');

%% 6. Thirs order intermodulation products of the amplifier

% Plot the 3rd order power

p_in = -20:2:0;         % Input power (dBm)
% Measured power of the 3rd intermodulation product
p_3rd = [-85 -80 -73 -68 -61.4 -54 -45 -34.8 -24.8 -15.7 -11];

figure('Color',[1 1 1]);

plot(p_in, p_3rd); hold on;
xlabel('P_{in} (dBm)');
xlim([p_in(1) p_in(end)]);
ylabel('P_{out, 3^{rd}} (dBm)');

set(gcf,'position',[100,100,400,300]);
set(gca,'FontSize',11, 'fontname','Times New Roman');

save_path = '../../Lab3_ActiveMWDevices/Images/';
saveas(gca, [save_path, 'p_out_3rd'],'epsc');

% Plot the intermodulation distance

% The values of the main tones were not measured, but it is known from
% previous exercises the values of the gain for different input powers, so
% it is taken the value of the gain, distributed equally on both tones.
% This can be done without considering p_3rd because p_tone >> p_third
G_lineal = s21_mod(1:800:end);
G_lineal = G_lineal(1:11);
p_tone = p_in + G_lineal.' - 3;
im3 = p_tone - p_3rd;   % Intermodulation distance

figure('Color',[1 1 1]);

plot(p_in, im3); hold on;
xlabel('P_{in} (dBm)');
xlim([p_in(1) p_in(end)]);
ylabel('IM_3 (dB)');

set(gcf,'position',[100,100,400,300]);
set(gca,'FontSize',11, 'fontname','Times New Roman');

save_path = '../../Lab3_ActiveMWDevices/Images/';
saveas(gca, [save_path, 'im3'],'epsc');

%% Calculate the third order intermodulation point
% The toi is desfined as the point in which the hypothetical third order
% matches the first order. It is needed a new P_in axis for the tendency
% line. This is going to be more accurate because is going to be used to
% detect the toi
p_in_larger = -20:0.1:25;         % Input power (dBm)

first_order = p_in + G_lineal.';
hyp_1st_order = p_in + G_lineal(1);
poly_1st = polyfit(p_in, hyp_1st_order, 1);
hyp_1st_order = polyval(poly_1st, p_in_larger);

% The third_order is calculated as the linear polynomial that best fits the
% data collected
third_order = p_3rd + 3;
poly_3rd = polyfit(p_in(end-1:end), third_order(end-1:end), 1);
hyp_3rd_order = polyval(poly_3rd, p_in_larger);

% The toi is at the point where the hyp_3rd > hyp_1st
loc_toi = find(hyp_3rd_order >= hyp_1st_order, 1);
fprintf('TOI = %f dBm\n', hyp_3rd_order(loc_toi));

figure('Color',[1 1 1]);

plot(p_in, first_order); hold on;
plot(p_in_larger, hyp_1st_order, '--'); hold on;
plot(p_in, third_order); hold on;
plot(p_in_larger, hyp_3rd_order, '--'); hold on;
plot(p_in_larger(loc_toi), hyp_3rd_order(loc_toi), 'o'); hold on;

xlim([p_in_larger(1) p_in_larger(end)]);
xlabel('Input power (dBm)');
ylabel('Output power (dBm)');
legend('P_{1^{st}}', 'Hypothetical P_{1^{st}}', 'P_{3^{rd}}', ...
    'Hypothetical P_{3^{rd}}', 'TOI', 'location', 'southeast');

set(gcf,'position',[100,100,400,300]);
set(gca,'FontSize',11, 'fontname','Times New Roman');

saveas(gca, [save_path, 'toi'],'epsc');

%% Find M-IMR
% In order to find ACPR it is necessary to open all the signals. All the
% files are opened. The frequency axis is common. The data is stored in a
% matrix nPoints x 11 because 11 is the number of measurements taken
n_points = 10001;
files = dir('Measurements_allTones/Multitone_Pin_*');
data = zeros(n_points, length(files));

P_in = -20:2:0;
% For each one of the files it is looked for the input power (in the name
% of the file). This is used to store the file in the adequate position
for ii = 1:length(files)
    filename = files(ii).name;
    input_power = str2double(filename(15:17));
    column_data = find(P_in == input_power, 1);
    
    % Obtain the data. There are 45 header lines
    measured_results = readtable(filename, 'HeaderLines',45);
    measured_results = measured_results{:,:};

    data(:,column_data) = measured_results(:,2);
end

f = measured_results(:,1);
% The position of the peaks are separated 10 MHz: it is only necessary to
% find the power at each position and compare with the maximum power. The
% peaks are separated 10 MHz, which are 100 points
f_peaks = f(1:100:end);
data_peaks = data(1:100:end, :);
% The M_IMR is obtained subtracting the power of the maximum tone to the
% power of the given tone
M_IMR = max(data_peaks) - data_peaks;

% The M_IMR is plotted as a fnction of the frequency
figure('Color',[1 1 1]);

plot(f_peaks/1e9, M_IMR, '-+');
xlabel('Frequency (GHz)');
ylabel('dB');
ylim([10 90]);

% The legend is created iteratively
txt = cell(length(P_in),1);
for ii = 1:length(P_in)
   txt{ii}= sprintf('P_{in} = %i dBm', P_in(ii));
end
legend(txt, 'location', 'eastoutside')

set(gcf,'position',[100,100,900,300]);
set(gca,'FontSize',11, 'fontname','Times New Roman');

saveas(gca, [save_path, 'm_imr_f'],'epsc');

% It is very difficult to plot all the imr as a function of the input
% power. Therefore, they are plotted graphs_plot frequencies for each plot.
% The legend is also generated iteratively
graphs_plot = 5;
figure('Color',[1 1 1]);
count_graph = 1;
leg = cell(length(graphs_plot),1);

for ii = 1:length(f_peaks)
    if f_peaks(ii) >= 2.486e9 && f_peaks(ii) <= 2.514e9
        continue;
    end
    
    plot(P_in, M_IMR(ii,:).'); hold on;
    
    if mod(ii, graphs_plot) == 0
        leg{graphs_plot} = sprintf('f = %.3f MHz', f_peaks(ii)/1e9);
    else
        leg{mod(ii, graphs_plot)} = sprintf('f = %.3f MHz', f_peaks(ii)/1e9);
    end
    
    % If it is time to start a new graph it is necessary to incude the 
    % legend, save the previous graph and start a new figure
    if mod(ii, graphs_plot) == graphs_plot -1
        xlabel('P_{in} (dBm)');
        ylabel('dB');
        
        legend(leg, 'location', 'northeast');
        
        set(gcf,'position',[100,100,300,300]);
        set(gca,'FontSize',11, 'fontname','Times New Roman');
        
        saveas(gca, [save_path, 'm_imr_pin_', num2str(count_graph)],'epsc');
        count_graph = count_graph+1;
        
        figure('Color',[1 1 1]);
    end
end
%% Find the ACPR
% This is done integrating the peaks in data_peaks depending on their
% frequency
plot(f_peaks, data_peaks)
f_lb_lim = 2.486e9;     % First central frequency of the central part
f_ub_lim = 2.514e9;     % Last central frequency of the central part
f_lb_lim_pos = find(f_peaks == f_lb_lim, 1);    % Position of 1st freq
f_ub_lim_pos = find(f_peaks == f_ub_lim, 1);    % Position of last freq

data_peaks_lb = data_peaks(1:f_lb_lim_pos-1, :);
data_peaks_mt = data_peaks(f_lb_lim_pos:f_ub_lim_pos, :);
data_peaks_ub = data_peaks(f_ub_lim_pos+1:end, :);

% Integrate the power in each band
% log10 takes the log of all the numbers
% sum sums each column. We are interested in adding up each column
pow_lb = 10*log10(sum(10.^(data_peaks_lb/10))); % Power of lower band (dBm)
pow_mt = 10*log10(sum(10.^(data_peaks_mt/10))); % Power of multitones (dBm)
pow_ub = 10*log10(sum(10.^(data_peaks_ub/10))); % Power of upper band (dBm)
    % Power of the lower band and the side band
pow_sb = 10*log10(sum(10.^([data_peaks_lb; data_peaks_ub]/10)));

% Calculate the ACPR
ACPR_t = pow_mt - pow_sb;
ACPR_l = pow_mt - pow_lb;
ACPR_u = pow_mt - pow_ub;

% Plot the results
figure('Color',[1 1 1]);

plot(P_in, ACPR_t); hold on;
plot(P_in, ACPR_l); hold on;
plot(P_in, ACPR_u); hold on;

xlabel('P_{in} (dBm)');
ylabel('dB');

legend('ACPR_T', 'ACPR_L', 'ACPR_U', 'location', 'northeast');

set(gcf,'position',[100,100,400,300]);
set(gca,'FontSize',11, 'fontname','Times New Roman');

saveas(gca, [save_path, 'acpr'],'epsc');

%% Multitone + notch: NPR
% In order to find NPR it is necessary to open all the signals. All the
% files are opened. The frequency axis is common. The data is stored in a
% matrix nPoints x 11 because 11 is the number of measurements taken
n_points = 10001;
files = dir('Measurements_allTones-1/Multitone_*');
data = zeros(n_points, length(files));

P_in = -20:2:0;
% For each one of the files it is looked for the input power (in the name
% of the file). This is used to store the file in the adequate position
for ii = 1:length(files)
    filename = files(ii).name;
    input_power = str2double(filename(20:22));
    column_data = find(P_in == input_power, 1);
    
    % Obtain the data. There are 45 header lines
    measured_results = readtable(filename, 'HeaderLines',45);
    measured_results = measured_results{:,:};

    data(:,column_data) = measured_results(:,2);
end

f = measured_results(:,1);
% The position of the peaks are separated 10 MHz: it is only necessary to
% find the power at each position and compare with the maximum power. The
% peaks are separated 10 MHz, which are 100 points. 
f_peaks = f(1:100:end);
data_peaks = data(1:100:end, :);
% It is only interesting the central part (the multitones)
f_peaks_mt = f_peaks(f_lb_lim_pos:f_ub_lim_pos);
data_peaks_mt = data_peaks(f_lb_lim_pos:f_ub_lim_pos, :);

% It is necessary to compare the power of the suppressed tone with respect
% the total power of the other harmonics
f_sup = 2.497e9;
loc_f_sup = find(f_peaks_mt == f_sup, 1);

% It is calculated the power of the other tones (dBm) and the NPR
pow_mt = 10*log10...
    (sum(10.^(data_peaks_mt([1:(loc_f_sup-1) (loc_f_sup+1):end],:)/10)));
NPR = pow_mt - data_peaks_mt(loc_f_sup);

% Plot the results
figure('Color',[1 1 1]);

plot(P_in, NPR); hold on;

xlabel('P_{in} (dBm)');
ylabel('NPR (dB)');

set(gcf,'position',[100,100,400,300]);
set(gca,'FontSize',11, 'fontname','Times New Roman');

saveas(gca, [save_path, 'npr'],'epsc');

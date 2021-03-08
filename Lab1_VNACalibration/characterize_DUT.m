function [] = characterize_DUT(dut, a, b, c, alpha, beta, gamma, r22_rho22)
save_path = '../../Lab1_VNACalibration/Images/';

% Open the measurements
[s_par_measured, freq] = getSParametersFromFile_Quadrupole(dut);

% Transform the S-parameters to T-parameters
R_measured = s_param_to_t_param(s_par_measured);

% Load the calibration
% load('calibration_TRL.mat');

% Modify the T-parameters using the calibration 
R_dut = obtain_R_DUT(a, b, c, alpha, beta, gamma, r22_rho22, R_measured);

% Transform the T-parameters to S-parameters
S_dut = t_param_to_s_param(R_dut);

% Plot the S-parameters of the device (magnitude)
s11 = 20*log10(S_dut(1, 1, :));        s11 = s11(:);
s21 = 20*log10(S_dut(2, 1, :));        s21 = s21(:);
s12 = 20*log10(S_dut(1, 2, :));        s12 = s12(:);
s22 = 20*log10(S_dut(2, 2, :));        s22 = s22(:);

figure('Color',[1 1 1]);
set(gcf,'position',[100,100,400,300]);
set(gca,'FontSize',11, 'fontname','Times New Roman');

plot(freq/1e9, s11); hold on;
plot(freq/1e9, s12); hold on;
plot(freq/1e9, s21); hold on;
plot(freq/1e9, s22); hold on;
xlim([freq(1) freq(end)]/1e9);
xlabel('Frequency (GHz)');
ylabel('dB');
legend('s_{11}', 's_{12}', 's_{21}', 's_{22}', 'Location', 'southeast');
saveas(gca, [save_path, dut(1:end-4), '_mag'],'epsc');

% Plot the S-parameters of the device (phase)
phase_s11 = angle(S_dut(1, 1, :));    phase_s11 = phase_s11(:)*180/pi;
phase_s12 = angle(S_dut(1, 2, :));    phase_s12 = phase_s12(:)*180/pi;
phase_s21 = angle(S_dut(2, 1, :));    phase_s21 = phase_s21(:)*180/pi;
phase_s22 = angle(S_dut(2, 2, :));    phase_s22 = phase_s22(:)*180/pi;

figure('Color',[1 1 1]);
set(gcf,'position',[100,100,400,300]);
set(gca,'FontSize',11, 'fontname','Times New Roman');

plot(freq/1e9, phase_s11); hold on;
plot(freq/1e9, phase_s12); hold on;
plot(freq/1e9, phase_s21); hold on;
plot(freq/1e9, phase_s22); hold on;
xlabel('Frequency (GHz)');
xlim([freq(1) freq(end)]/1e9);
ylabel('Phase (º)');
legend('<s_{11}', '<s_{12}', '<s_{21}', '<s_{22}', 'Location', 'best');
saveas(gca, [save_path, dut(1:end-4), '_phase'],'epsc');
end


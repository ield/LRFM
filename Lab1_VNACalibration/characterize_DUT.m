% This script characterizes an object measured in the vna wthout
% calibration. It should be run after calibrationParameters.m so that the
% vna is calibrated correctly
clear;
close all;
save_path = '../../Lab1_VNACalibration/Images/';

% Open the measurements
[s_obj_measured] = getSParametersFromFile_Quadrupole('test_t.s2p');

% Transform the S-parameters to T-parameters
R_measured = s2t(s_obj_measured.Parameters);

% Load the calibration
load('calibration_TRL.mat');

% Modify the T-parameters using the calibration 
R_dut = obtain_R_DUT(a, b, c, alpha, beta, gamma, r22rho22, R_measured);

% Transform the T-parameters to S-parameters
S_dut = t2s(R_dut);

% Plot the S-parameters of the device (magnitude)
obj_s_dut = sparameters(S_dut, s_obj_measured.Frequencies, ...
    s_obj_measured.Impedance);
figure('Color',[1 1 1]);
set(gcf,'position',[100,100,400,300]);
set(gca,'FontSize',11, 'fontname','Times New Roman');
rfplot(obj_s_dut)
xlim([obj_s_dut.Frequencies(1) obj_s_dut.Frequencies(end)]/1e9);
saveas(gca, [save_path, 'dut_spar_mag'],'epsc');

% Plot the S-parameters of the device (phase)
phase_s11 = angle(S_dut(1, 1, :));    phase_s11 = phase_s11(:)*180/pi;
phase_s12 = angle(S_dut(1, 2, :));    phase_s12 = phase_s12(:)*180/pi;
phase_s21 = angle(S_dut(2, 1, :));    phase_s21 = phase_s21(:)*180/pi;
phase_s22 = angle(S_dut(2, 2, :));    phase_s22 = phase_s22(:)*180/pi;

figure('Color',[1 1 1]);
set(gcf,'position',[100,100,400,300]);
set(gca,'FontSize',11, 'fontname','Times New Roman');

plot(obj_s_dut.Frequencies/1e9, phase_s11); hold on;
plot(obj_s_dut.Frequencies/1e9, phase_s12); hold on;
plot(obj_s_dut.Frequencies/1e9, phase_s21); hold on;
plot(obj_s_dut.Frequencies/1e9, phase_s22); hold on;
xlim([obj_s_dut.Frequencies(1) obj_s_dut.Frequencies(end)]/1e9);
legend('<s_{11}', '<s_{12}', '<s_{21}', '<s_{22}', 'Location', 'best');
saveas(gca, [save_path, 'dut_spar_phase'],'epsc');
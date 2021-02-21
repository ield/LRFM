% This script characterizes an object measured in the vna wthout
% calibration. It should be run after calibrationParameters.m so that the
% vna is calibrated correctly
clear;
close all;

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

% Plot the S-parameters of the device
obj_s_dut = sparameters(S_dut, s_obj_measured.Frequencies, ...
    s_obj_measured.Impedance);
figure('Color',[1 1 1]);
set(gcf,'position',[100,100,400,300]);
set(gca,'FontSize',11, 'fontname','Times New Roman');
rfplot(obj_s_dut)
xlim([obj_s_dut.Frequencies(1) obj_s_dut.Frequencies(end)]/1e9);
title('S-parameters of DUT');
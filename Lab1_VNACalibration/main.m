% Complete the practice in MATLAB

% Calibrate using line 1
calibrate_TRL('thru.dat', 'line_1_135mm.dat', 'reflect_gate1.dat', ...
    'reflect_gate2.dat');
load('calibration_TRL.mat');

%%

characterize_DUT('dut.dat', a, b, c, alpha, beta, gamma, r22_rho22);
characterize_DUT('line_1_135mm.dat', a, b, c, alpha, beta, gamma, ...
    r22_rho22);
characterize_DUT('line_2_270mm.dat', a, b, c, alpha, beta, gamma, ...
    r22_rho22);

%%
obtain_delay('line_1_135mm.dat', a, b, c, alpha, beta, gamma, r22_rho22);
obtain_delay('line_2_270mm.dat', a, b, c, alpha, beta, gamma, r22_rho22);

%%
characterize_reflect('reflect_gate1.dat', 'reflect_gate2.dat', ...
    a, b, c, alpha, beta, gamma, freq);

%% Calibrate using line2
clear;
calibrate_TRL('thru.dat', 'line_2_270mm.dat', 'reflect_gate1.dat', ...
    'reflect_gate2.dat');    
load('calibration_TRL.mat');
%% 
characterize_DUT('line_1_135mm.dat', a, b, c, alpha, beta, gamma, ...
    r22_rho22);
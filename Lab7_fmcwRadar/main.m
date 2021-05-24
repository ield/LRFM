clc;clear; close all;

%% General data
c = 2e8;
f_m = 135.63;
Delta_f_1 = 120e6;
Delta_f_2 = 60e6;
%% Exercise 1
% Time domain
% Datos 1
x1 = 9.216e-4;
x2 = 5.184e-3;
n_periods = 6;

T_b = (x2-x1)/n_periods;
f_b = 1/T_b;

R = c*f_b/(f_m*Delta_f_1);

fprintf('Time Domain, datos1. T_b = %f ms. f_b = %f kHz. R = %f m\n', ...
    T_b*1e3, f_b/1e3, R);

% Datos 2
x1 = 1.728e-3;
x2 = 5.9328e-3;
n_periods = 3;

T_b = (x2-x1)/n_periods;
f_b = 1/T_b;

R = c*f_b/(f_m*Delta_f_2);

fprintf('Time Domain, datos2. T_b = %f ms. f_b = %f kHz. R = %f m\n', ...
    T_b*1e3, f_b/1e3, R);

% Frequency domain
% Datos 1
f_b = 1.4215e3;

R = c*f_b/(f_m*Delta_f_1);

fprintf('Frequency Domain, datos1. R = %f m\n', R);

% Datos 2
f_b = 0.678168e3;

R = c*f_b/(f_m*Delta_f_2);

fprintf('Frequency Domain, datos2. R = %f m\n', R);

%% Sampling frequency
f_b_max = 8.33e3;
f_s = 2*f_b_max;
fprintf('f_s = %f kHz\n', f_s/1e3);

R_max = c*f_b_max/(f_m*Delta_f_1);

fprintf('datos1. Rmax = %f m\n', R_max);

R_max = c*f_b_max/(f_m*Delta_f_2);

fprintf('datos2. Rmax = %f m\n', R_max);

%% Noise bandwidth
s_on = -21
s_off = -20;

diff = s_off - s_on;
b_noise = 10^((10*log10(f_b_max) + diff)/10);
fprintf('B = %f kHz\n', b_noise/1e3);

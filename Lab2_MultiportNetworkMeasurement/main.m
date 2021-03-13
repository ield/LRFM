% Author: ield
% Subject: LRFM MSTC ETSIT UPM
%% Description of the code
% This code gathers the work done in order to characterize a multiport
% network.
% 1. All the S-parameters are taken from the measurements done in the Lab
% 2. The S-parameters measured are rerreferred to s-parameters of the
% device
% 3. The S-parameters are transformed from Z0 to Zl
% 4. The repeated values are checked to be equal.
% 5. The S-parameters are transformed from Z0 to Zl
clear;
clc;
close all;
save_path = '../../Lab2_MultiportNetworkMeasurement/Images/';

%% Import S-parameters
% The naming is ae[port of of the ae] __connected to__ dut[port of the dut]
% l[port of the load]
[s_par_ae1_dut1_ae2_dut2_l3, f] = ...
    getSParametersFromFile_Quadrupole('ae1_dut1_ae2_dut2_l3.dat');
[s_par_ae1_dut1_ae2_dut3_l2, ~] = ...
    getSParametersFromFile_Quadrupole('ae1_dut1_ae2_dut3_l2.dat');
[s_par_ae1_dut2_ae2_dut3_l1, ~] = ...
    getSParametersFromFile_Quadrupole('ae1_dut2_ae2_dut3_l1.dat');
[s11_load, ~] = getSParametersFromFile_Reflection('load_thru.dat');

% The reference impedance is 30 Ohm
Z0 = 30;

% Plot the S-parameters of the measurement.
plotS_Parameters_2port(s_par_ae1_dut1_ae2_dut2_l3, f, save_path, ...
    'meas_1_z0');
plotS_Parameters_2port(s_par_ae1_dut1_ae2_dut3_l2, f, save_path, ...
    'meas_2_z0');
plotS_Parameters_2port(s_par_ae1_dut2_ae2_dut3_l1, f, save_path, ...
    'meas_3_z0');
%% Modify the reference impedance of the S-parameters
% The s-parameters are transformed from z0 tozl1, zl2, zl3... Since zl1 =
% zl2 = zl3, they are all transformed to zl

% In order to transform from Z0 to Zl first it is necessary to characterize
% Zl
Zl = Z0*(1+s11_load)./(1-s11_load);

% Then it s followed the proceduredescribed in slide 31, with Zb = Zl and
% Za = Z0. It is necessary to define the new sparameters with respect to zl
s_par_ae1_dut1_ae2_dut2_l3_zl = zeros(size(s_par_ae1_dut1_ae2_dut2_l3));
s_par_ae1_dut1_ae2_dut3_l2_zl = zeros(size(s_par_ae1_dut1_ae2_dut3_l2));
s_par_ae1_dut2_ae2_dut3_l1_zl = zeros(size(s_par_ae1_dut2_ae2_dut3_l1));

% The operation mus be carried out for all f
for ii = 1:length(f)
    eta = (Zl(ii) - Z0)/(Z0+Zl(ii))*eye(2);
    G = sqrt(Z0/Zl(ii))*eye(2);
    U = eye(2);
    
    % For the measure s_par_ae1_dut1_ae2_dut2_l3
    Sa = s_par_ae1_dut1_ae2_dut2_l3(:,:,ii);
    s_par_ae1_dut1_ae2_dut2_l3_zl(:,:,ii) = ...
        G*inv(U-Sa)*(Sa-eta)*inv(U-Sa*eta)*(U-Sa)*inv(G);
    
    % For the measure s_par_ae1_dut1_ae2_dut3_l2
    Sa = s_par_ae1_dut1_ae2_dut3_l2(:,:,ii);
    s_par_ae1_dut1_ae2_dut3_l2_zl(:,:,ii) = ...
        G*inv(U-Sa)*(Sa-eta)*inv(U-Sa*eta)*(U-Sa)*inv(G);
    
    % For the measure s_par_ae1_dut2_ae2_dut3_l1
    Sa = s_par_ae1_dut2_ae2_dut3_l1(:,:,ii);
    s_par_ae1_dut2_ae2_dut3_l1_zl(:,:,ii) = ...
        G*inv(U-Sa)*(Sa-eta)*inv(U-Sa*eta)*(U-Sa)*inv(G);
end

% Plot the new S-parameters referred to Zl
plotS_Parameters_2port(s_par_ae1_dut1_ae2_dut2_l3_zl, f, save_path, ...
    'meas_1_zl');
plotS_Parameters_2port(s_par_ae1_dut1_ae2_dut3_l2_zl, f, save_path, ...
    'meas_2_zl');
plotS_Parameters_2port(s_par_ae1_dut2_ae2_dut3_l1_zl, f, save_path, ...
    'meas_3_zl');
%% Refer the S-parameters to the device
% There are some of the parameters duplicated. They are differentiated by
% the load port.
% After they are taken, they are converted to single vector form

%   s_par_ae1_dut1_ae2_dut2_l3
%       s11 is the s11 of the device. This will be duplicated
%       s12 is the s12 of the device
%       s21 is the s21 of the device
%       s22 is the s22 of the device. This will be duplicated
s11_l3_zl = s_par_ae1_dut1_ae2_dut2_l3_zl(1, 1, :);    
s11_l3_zl = s11_l3_zl(:);
s12_zl = s_par_ae1_dut1_ae2_dut2_l3_zl(1, 2, :);       
s12_zl = s12_zl(:);
s21_zl = s_par_ae1_dut1_ae2_dut2_l3_zl(2, 1, :);       
s21_zl = s21_zl(:);
s22_l3_zl = s_par_ae1_dut1_ae2_dut2_l3_zl(2, 2, :);    
s22_l3_zl = s22_l3_zl(:);

%   s_par_ae1_dut1_ae2_dut3_l2
%       s11 is the s11 of the device. This will be duplicated
%       s12 is the s13 of the device
%       s21 is the s31 of the device
%       s22 is the s33 of the device. This will be duplicated
s11_l2_zl = s_par_ae1_dut1_ae2_dut3_l2_zl(1, 1, :);    
s11_l2_zl = s11_l2_zl(:);
s13_zl = s_par_ae1_dut1_ae2_dut3_l2_zl(1, 2, :);       
s13_zl = s13_zl(:);
s31_zl = s_par_ae1_dut1_ae2_dut3_l2_zl(2, 1, :);       
s31_zl = s31_zl(:);
s33_l2_zl = s_par_ae1_dut1_ae2_dut3_l2_zl(2, 2, :);    
s33_l2_zl = s33_l2_zl(:);

%   s_par_ae1_dut2_ae2_dut3_l1
%       s11 is the s22 of the device. This will be duplicated
%       s12 is the s23 of the device
%       s21 is the s32 of the device
%       s22 is the s33 of the device. This will be duplicated
s22_l1_zl = s_par_ae1_dut2_ae2_dut3_l1_zl(1, 1, :);    
s22_l1_zl = s22_l1_zl(:);
s23_zl = s_par_ae1_dut2_ae2_dut3_l1_zl(1, 2, :);       
s23_zl = s23_zl(:);
s32_zl = s_par_ae1_dut2_ae2_dut3_l1_zl(2, 1, :);       
s32_zl = s32_zl(:);
s33_l1_zl = s_par_ae1_dut2_ae2_dut3_l1_zl(2, 2, :);    
s33_l1_zl = s33_l1_zl(:);

% View the error between the redundant parameters
redundant_parameters(s11_l3_zl, s11_l2_zl, f, 'Meas 1', 'Meas 2', ...
    save_path, 's11_redundant');
redundant_parameters(s22_l3_zl, s22_l1_zl, f, 'Meas 1', 'Meas 3', ...
    save_path, 's22_redundant');
redundant_parameters(s33_l2_zl, s33_l1_zl, f, 'Meas 2', 'Meas 3', ...
    save_path, 's33_redundant');

%% Transform S-parameters from Z0 to Zl
s_par_zl = zeros(3, 3, length(f));
s_par_zl(1, 1, :) = s11_l3_zl;
s_par_zl(1, 2, :) = s12_zl;
s_par_zl(1, 3, :) = s13_zl;
s_par_zl(2, 1, :) = s21_zl;
s_par_zl(2, 2, :) = s22_l3_zl;
s_par_zl(2, 3, :) = s23_zl;
s_par_zl(3, 1, :) = s31_zl;
s_par_zl(3, 2, :) = s32_zl;
s_par_zl(3, 3, :) = s33_l1_zl;

% Once all the parameters are taken, they are plotted in mod and phase
plotS_Parameters_3port(s_par_zl, f, save_path, 's_par_zl');

s_par_z0 = zeros(3, 3, length(f));

for ii = 1:length(f)
    eta = (Z0 - Zl(ii))/(Z0 + Zl(ii))*eye(3);
    G = sqrt(Zl(ii)/Z0)*eye(3);
    U = eye(3);
    Sa = s_par_zl(:,:,ii);
    
    s_par_z0(:,:,ii) = G*inv(U-Sa)*(Sa-eta)*inv(U-Sa*eta)*(U-Sa)*inv(G);
end

% Once all the parameters are taken, they are plotted in mod and phase
plotS_Parameters_3port(s_par_z0, f, save_path, 's_par_z0');

%% Figure out the purpose of the device
% It is checked for symmetries: s_mn = s_nm
parToPlot = [1 2; 2 1; 1 3; 3 1; 2 3; 3 2];
plot_Transmission_Parameters(s_par_z0, parToPlot, f, save_path, ...
    'reciprocal');

% It is checked if the device is a -3 dB hybrid coupler
parToPlot = [1 1; 2 1; 3 1];
plot_Transmission_Parameters(s_par_z0, parToPlot, f, save_path, ...
    'hybrid', 1);




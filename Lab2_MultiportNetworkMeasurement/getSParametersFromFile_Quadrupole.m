%Engineer: ield

% function [f, s11, s21, s12, s22] = getSParametersFromFile_Quadrupole(...
%     filename)
function [s_par, f] = getSParametersFromFile_Quadrupole(filename)
%% General Information
% Reads the .dat file generated with the VNA and returns the f, and the
% s-parameters

% It reads the files
file_info = importdata(filename);

freq_unit = 1;      % Units of frequency (respect to Hz). Hz = 1. GHz = 1e9
Z0 = 50;            % Reference impedance, Ohm

f = file_info.data(:,1)/freq_unit;

% The order of the s parameters is s21, s11, s12, s22 real and imag
s_par = zeros(2, 2, length(f));
s_par(1, 1, :) = file_info.data(:,4) + 1j*file_info.data(:,5);
s_par(2, 1, :) = file_info.data(:,2) + 1j*file_info.data(:,3);
s_par(1, 2, :) = file_info.data(:,6) + 1j*file_info.data(:,7);
s_par(2, 2, :) = file_info.data(:,8) + 1j*file_info.data(:,9);

end


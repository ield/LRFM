%Engineer: ield

function [s_11, f] = getSParametersFromFile_Reflection(filename)
%% General Information
% Reads the .dat file generated with the VNA and returns the f, and the
% s-parameters

% It reads the files
file_info = importdata(filename);

freq_unit = 1;      % Units of frequency (respect to Hz). Hz = 1. GHz = 1e9
Z0 = 50;            % Reference impedance, Ohm

f = file_info.data(:,1)/freq_unit;

% s_11 = zeros(1, length(f));
s_11 = file_info.data(:,2) + 1j*file_info.data(:,3);
end


%Engineer: ield

function [s_obj] = getSParametersFromFile_Reflection(filename)
%% General Information
% Reads the .s2p file generated with the VNA and returns the f, and the
% s-parameters

% It reads the files
file_spar = fopen(filename);
data_spar = textscan(file_spar,'%s');
fclose(file_spar);
initial_data = 16;  % Line of the first frequency evaluated
data_spar = str2double(data_spar{1}(initial_data:end));

freq_unit = 1;      % Units of frequency (respect to Hz). Hz = 1. GHz = 1e9
Z0 = 50;            % Reference impedance, Ohm
num_columns = 3;    % Columns: f + real and imag of each sparameter

f = data_spar(1:num_columns:end)/freq_unit;

s_par = zeros(1, 1, length(f));
s_par(1, 1, :) = data_spar(2:num_columns:end) + 1j*data_spar(3:num_columns:end);

s_obj = sparameters(s_par, f, Z0);
end


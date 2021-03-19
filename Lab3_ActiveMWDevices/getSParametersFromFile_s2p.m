%Engineer: ield

function [axis, s11_mod, s11_phase, s21_mod, s21_phase, s12_mod, s12_phase, ...
    s22_mod, s22_phase] = getSParametersFromFile_s2p(filename)
%% General Information
% Reads the .s2p file generated with the VNA and returns the f, and the
% s-parameters

% It reads the files
fmt = repmat('%f', 1, 9);
fid = fopen(filename, 'rt');
datacell = textscan(fid, fmt, 'HeaderLines', 6, 'CollectOutput', 1);
fclose(fid);
data = datacell{1};

axis = data(:,1);
s11_mod = data(:,2);
s11_phase = data(:,3);
s21_mod = data(:,4);
s21_phase = data(:,5);
s12_mod = data(:,6);
s12_phase = data(:,7);
s22_mod = data(:,8);
s22_phase = data(:,9);

end


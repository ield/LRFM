function [data_filtered] = lp_filter_data(data)
% Low-pass filtering of a given data
spectrum_data = fft(data);

% It is created a pseudo frequency domain
L = length(data);
f = (0:(L/2))/L;

f_cut = 1/L;
index_f_cut = find (f>f_cut, 1);

% It is created a ideal low pass filter of fc= f-cut_
lp_filter = zeros(size(spectrum_data));
lp_filter(1:index_f_cut) = ones(size(lp_filter(1:index_f_cut)));
lp_filter(end-(index_f_cut-2):end) = ...
    ones(size(lp_filter(end-(index_f_cut-2):end)));

% The signal is filtered
spectrum_data = spectrum_data.*lp_filter;
data_filtered = ifft(spectrum_data);
end


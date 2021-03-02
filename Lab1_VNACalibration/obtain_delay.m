function [] = obtain_delay(line, a, b, c, alpha, beta, gamma, r22_rho22)
save_path = '../../Lab1_VNACalibration/Images/';

% Open the measurements
[s_line, freq] = getSParametersFromFile_Quadrupole(line);

% Transform the S-parameters to T-parameters
R_L = s_param_to_t_param(s_line);

% Obtain the S-parameters of the line
R_line = obtain_R_DUT(a, b, c, alpha, beta, gamma, r22_rho22, R_L);
S_line = t_param_to_s_param(R_line);

% Obtain the phase delay of the line
s21_line = S_line(2,1,:);   s21_line = s21_line(:);
length_line = str2double(line(8:10))*1e-4;      % Length of the line (m)
fprintf('The length of the line is %i mm\n', length_line*1000);
gamma_m = log(s21_line)/(-length_line);
beta_m = imag(gamma_m);
phase_delay = beta_m*length_line./(2*pi*freq);

% Plot the phase delay of the line
figure('Color',[1 1 1]);
set(gcf,'position',[100,100,400,300]);
set(gca,'FontSize',11, 'fontname','Times New Roman');
plot(freq/1e9, phase_delay*1e12);
title('Phase delay of the line(\tau_l)');
xlabel('f (GHz)');
ylabel('ps');
xlim([freq(1) freq(end)]/1e9);
saveas(gca, [save_path, line(1:end-4), '_delay'],'epsc');

fprintf('The mean phase delay is %f s\n', mean(phase_delay(2:300)*1e12));
fmin = 3e8 / (18*length_line);
fmax = 4*3e8 / (9*length_line);

fprintf('fmin = %f GHz, fmax = %f GHz\n', fmin/1e9, fmax/1e9);

end


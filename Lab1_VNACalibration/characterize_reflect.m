function [] = characterize_reflect(ref_1, ref_2, ...
    a, b, c, alpha, beta, gamma, freq)

save_path = '../../Lab1_VNACalibration/Images/';

% Open the measurements of the reflections
[r1_s11] = getSParametersFromFile_Reflection(ref_1);    r1_s11 = r1_s11(:);
[r2_s11] = getSParametersFromFile_Reflection(ref_2);    r2_s11 = r2_s11(:);

% Calculate both reflection coefficientes of the reflect and see the
% differences
rho_R_1 = (r1_s11 - b)./(a - r1_s11.*c);
rho_R_2 = (gamma + r2_s11)./(alpha + beta.*r2_s11);
error = abs(rho_R_1 - rho_R_2);

% Plot reflections and the error
figure('Color',[1 1 1]);
set(gcf,'position',[100,100,1200,300]);
set(gca,'FontSize',11, 'fontname','Times New Roman');

subplot(1, 3, 1)
smithplot(freq, rho_R_1, 'TitleTop', 'Reflection port 1');
title('Reflect on gate 1');

subplot(1, 3, 2)
smithplot(freq, rho_R_2, 'TitleTop', 'Reflection port 2');
title('Reflect on gate 2');

subplot(1, 3, 3);
plot(freq, error);
xlim([freq(1) freq(end)]/1e9);
title('|\Gamma_{R port 1} - \Gamma_{R port 2}|');
xlabel('f (GHz)');
% saveas(gca, [save_path, 'reflection'],'epsc');

% Plot the magnitude and phase of the reflect
figure('Color',[1 1 1]);
set(gcf,'position',[100,100,1200,300]);
set(gca,'FontSize',11, 'fontname','Times New Roman');

subplot(1, 2, 1)
plot(freq/1e9, -20*log10(abs(rho_R_1)));    hold on;
plot(freq/1e9, -20*log10(abs(rho_R_2)));    
xlim([freq(1) freq(end)]/1e9);
xlabel('f (GHz)');
ylabel('dB');
title('Magnitude');
legend('Gate 1', 'Gate 2');

subplot(1, 2, 2)
plot(freq/1e9, angle(rho_R_1)*180/pi);    hold on;
plot(freq/1e9, angle(rho_R_2)*180/pi);    hold on;
xlim([freq(1) freq(end)]/1e9);
title('Phase');
xlabel('f (GHz)');
ylabel('(º)');
legend('Gate 1', 'Gate 2');

saveas(gca, [save_path, 'reflection_mag_phase'],'epsc');
end


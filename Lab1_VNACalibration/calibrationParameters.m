clear;
save_path = '../../Lab1_VNACalibration/Images/';
%% Open the measurements
% Open the S-Parameters of the thru and the line, measures taken with the
% VNA
[s_obj_t] = getSParametersFromFile_Quadrupole('test_t.s2p');
[s_obj_l] = getSParametersFromFile_Quadrupole('test_l.s2p');

% Open the measurements of the reflections
[s_obj_r1] = getSParametersFromFile_Reflection('test_r.s2p');
[s_obj_r2] = getSParametersFromFile_Reflection('test_r.s2p');

%% Steps on slide 20
% Transform the S-parameters to T-parameters (using the names in slide20)
R_T = s2t(s_obj_t.Parameters);
R_L = s2t(s_obj_l.Parameters);

% Calculate T
T = zeros(size(R_T));
for ii=1:size(T,3)
   T(:,:,ii) = R_L(:,:,ii) * inv(R_T(:,:,ii)); 
end

%% Steps on slide 21
% It is solved the second degree equation using the traditional formula
% -b+-sqrt... where a = t21, b = t22-t11, c = -t12
a_par = T(2,1,:);
b_par = T(2,2,:) - T(1,1,:);
c_par = -T(1,2,:);
sol_1 = (-b_par+sqrt(b_par.^2-4*a_par.*c_par))./(2.*a_par);
sol_2 = (-b_par-sqrt(b_par.^2-4*a_par.*c_par))./(2.*a_par);

% Now it is determines between the two solutions which one fits better the
% specifications: abs(b) << 1 and abs(a/c) >> 1. This is going to be done
% by calculatig the mean of the abslute values
if mean(abs(sol_1)) > mean(abs(sol_2))
    b = sol_2;
    a_over_c = sol_1;
else
    b = sol_1;
    a_over_c = sol_2;
end
fprintf('The mean of abs(b) is %f anf the mean of abs(a/c) is %f\n', ...
    mean(abs(b)), mean(abs(a_over_c)));

%% Steps on slide 22
% d, e, f, and g are extracted from the R_T matrix
g = R_T(2,2,:);
d = R_T(1,1,:)./g;
e = R_T(1,2,:)./g;
f = R_T(2,1,:)./g;

c_over_a = 1./a_over_c;

r22_rho22 = g.*(1-e.*c_over_a)./(1-b.*c_over_a);
gamma = (f-d.*c_over_a)./(1-e.*c_over_a);
alpha_a = (d - b.*f)./(1-e.*c_over_a);
beta_over_alpha = (e-b)./(d-b.*f);

%% Steps on slide 24
a_over_alpha = (s_obj_r1.Parameters(1,1,:) - b).*...
    (1 + s_obj_r2.Parameters(1,1,:).*beta_over_alpha)./...
    ((s_obj_r2.Parameters(1,1,:) + gamma).*...
    (1-s_obj_r1.Parameters(1,1,:).*c_over_a));

% There are two possible values of a. Both are calculated and plotted, then
% they are checked.
a_sol1 = sqrt(alpha_a .* a_over_alpha);
a_sol2 = -sqrt(alpha_a .* a_over_alpha);

rho_R_a_sol1 = (s_obj_r1.Parameters(1,1,:) - b)./...
    (a_sol1.*(1-s_obj_r1.Parameters(1,1,:).*c_over_a));
rho_R_a_sol2 = (s_obj_r1.Parameters(1,1,:) - b)./...
    (a_sol2.*(1-s_obj_r1.Parameters(1,1,:).*c_over_a));

% Plot both parameters
figure('Color',[1 1 1]);
set(gcf,'position',[100,100,400,300]);
set(gca,'FontSize',11, 'fontname','Times New Roman');
smithplot(s_obj_r1.Frequencies, rho_R_a_sol1(:)); hold on;
smithplot(s_obj_r1.Frequencies, rho_R_a_sol2(:)); hold on;
title('\Gamma_R for each solution of a');
legend('+ sol', '- sol', 'Location', 'eastoutside');
% saveas(gca, [save_path, 'smit_refflection'],'epsc');

% Check the best parameter: the one closest to one, because the reflection
% is an open circuit
error_rho_R_a_sol1 = sum(abs(rho_R_a_sol1 - 1).^2);
error_rho_R_a_sol2 = sum(abs(rho_R_a_sol2 - 1).^2);

if error_rho_R_a_sol1 < error_rho_R_a_sol2
    a = a_sol1;
    fprintf('It is taken the positive answer\n');
else
    a = a_sol2;
    fprintf('It is taken the negative answer\n');
end

% Now it is calculated c, alpha and beta
c = a.*c_over_a;
alpha = alpha_a ./ a;
beta = alpha .* beta_over_alpha;

%% Steps on slide 25
% Calculate both reflection coefficientes of the reflect and see the
% differences
rho_R_1 = (s_obj_r1.Parameters(1,1,:) - b)./...
    (a - s_obj_r1.Parameters(1,1,:).*c);
rho_R_2 = (gamma + s_obj_r2.Parameters(1,1,:))./...
    (alpha + beta.*s_obj_r2.Parameters(1,1,:));
error = abs(rho_R_1 - rho_R_2);

% Plot the error
figure('Color',[1 1 1]);
set(gcf,'position',[100,100,400,300]);
set(gca,'FontSize',11, 'fontname','Times New Roman');
plot(s_obj_r1.Frequencies/1e9, error(:));
xlim([s_obj_r1.Frequencies(1) s_obj_r1.Frequencies(end)]/1e9);
% title('|\Gamma_{R port 1} - \Gamma_{R port 2}|');
xlabel('f (GHz)');
saveas(gca, [save_path, 'error_reflection'],'epsc');

% Obtain the S-parameters of the line
R_line = obtain_R_DUT(a, b, c, alpha, beta, gamma, r22_rho22, R_L);
S_line = t2s(R_line);

% Obtain the phase delay of the line
s21_line = S_line(2,1,:);   s21_line = s21_line(:);
length_line = 10e-3;      % Length of the line (m)
gamma_m = log(s21_line)/(-length_line);
beta_m = imag(gamma_m);
phase_delay = beta_m*length_line./(2*pi*s_obj_r1.Frequencies);

% Plot the phase delay of the line
figure('Color',[1 1 1]);
set(gcf,'position',[100,100,400,300]);
set(gca,'FontSize',11, 'fontname','Times New Roman');
plot(s_obj_r1.Frequencies/1e9, phase_delay);
title('Phase delay of the line(\tau_l)');
xlabel('f (GHz)');
ylabel('s');
xlim([s_obj_r1.Frequencies(1) s_obj_r1.Frequencies(end)]/1e9);
saveas(gca, [save_path, 'phase_delay'],'epsc');

fprintf('The mean phase delay is %f s\n', mean(phase_delay));

%% Check that all the parameters fit the desired expectations
figure('Color',[1 1 1]);
set(gcf,'position',[100,100,400,300]);
set(gca,'FontSize',11, 'fontname','Times New Roman');
plot(s_obj_r1.Frequencies/1e9, abs(a(:))); hold on;
plot(s_obj_r1.Frequencies/1e9, abs(b(:))); hold on; 
plot(s_obj_r1.Frequencies/1e9, abs(c(:))); hold on;
plot(s_obj_r1.Frequencies/1e9, abs(alpha(:))); hold on;
plot(s_obj_r1.Frequencies/1e9, abs(beta(:))); hold on;
plot(s_obj_r1.Frequencies/1e9, abs(gamma(:))); hold on;
plot(s_obj_r1.Frequencies/1e9, abs(r22_rho22(:)));
xlabel('f (GHz)');
xlim([s_obj_r1.Frequencies(1) s_obj_r1.Frequencies(end)]/1e9);
legend('|a|', '|b|', '|c|', '|\alpha|', '|\beta|', '|\gamma|', ...
    '|r_{22} \rho_{22}|', 'Location', 'best');
saveas(gca, [save_path, 'abs_parameters'],'epsc');
%% Save all the calibration parameters
freq = s_obj_r1.Frequencies;
filename = 'calibration_TRL.mat';
save(filename, 'a', 'b', 'c', 'alpha', 'beta', 'gamma', 'r22rho22', ...
    'length_line', 'phase_delay', 'freq');

    
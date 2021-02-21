function [R_DUT] = obtain_R_DUT(a, b, c, alpha, beta, gamma, r22rho22, R_M)
%The function obtains the T matrix of the device under test with the
%parameters obtained from the TRL calibration and the T parameters of the
%measured. It is used the last expression in slide 19
R_DUT = zeros(size(R_M));

% Now they are created the matrices left and right of rm
left_rm = zeros(size(R_M));
left_rm(1, 1, :) = ones(size(R_M, 3), 1);
left_rm(1, 2, :) = -b(:);
left_rm(2, 1, :) = -c(:);
left_rm(2, 2, :) = a(:);

right_rm = zeros(size(R_M));
right_rm(1, 1, :) = ones(size(R_M, 3), 1);
right_rm(1, 2, :) = -beta(:);
right_rm(2, 1, :) = -gamma(:);
right_rm(2, 2, :) = alpha(:);

% Now it is completed the expression in slide 19. First it is calculated
% the first term
first_term = 1./(r22rho22.*(a-b.*c).*(alpha-beta.*gamma));
first_term = first_term(:);     % It is converted to a column
for ii = 1:size(R_M, 3)
    R_DUT(:,:,ii) = first_term(ii)*left_rm(:,:,ii)*R_M(:,:,ii)*...
        right_rm(:,:,ii);
end

end


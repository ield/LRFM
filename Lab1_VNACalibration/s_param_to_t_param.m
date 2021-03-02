function [t_param] = s_param_to_t_param(s_param)
    % Trasnforms from s-params to t-params according to the slides of the
    % subject lrfm. It is better than s2t matlab's function
    t_param = zeros(size(s_param));
    
    delta_s = s_param(1, 1, :).*s_param(2, 2, :) - ....
        s_param(2, 1, :).*s_param(1, 2, :);
    
    t_param(1, 1, :) = -delta_s./s_param(2, 1, :);
    t_param(1, 2, :) = s_param(1, 1, :)./s_param(2, 1, :);
    t_param(2, 1, :) = -s_param(2, 2, :)./s_param(2, 1, :);
    t_param(2, 2, :) = 1./s_param(2, 1, :);

end


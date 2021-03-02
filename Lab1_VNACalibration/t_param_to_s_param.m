function [s_param] = t_param_to_s_param(t_param)
    s_param = zeros(size(t_param));
    
    delta_t = t_param(1, 1, :).*t_param(2, 2, :) - ...
        t_param(2, 1, :).*t_param(1, 2, :);
    
    s_param(1, 1, :) = t_param(1, 2, :)./t_param(2, 2, :);
    s_param(1, 2, :) = delta_t./t_param(2, 2, :);
    s_param(2, 1, :) = 1./t_param(2, 2, :);
    s_param(2, 2, :) = -t_param(2, 1, :)./t_param(2, 2, :);

end


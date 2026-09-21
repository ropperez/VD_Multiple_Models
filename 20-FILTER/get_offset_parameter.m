function [b] = get_offset_parameter(r_k,r_k_1,a_cte,w_cte,ENUM)
% b: control parameter

% State transition matrix
enum_dyn_model = ENUM.dynamic_model;

% Mode in r_k_1
if r_k == enum_dyn_model.cv
    b   = [0; 0; 0; 0];

elseif r_k == enum_dyn_model.ca
    % Mode in r_k
    if r_k_1 == enum_dyn_model.cv
        b = [0; 0; a_cte.x; 0; 0; a_cte.y];

    elseif r_k_1 == enum_dyn_model.ca
        b = [0; 0; 0; 0; 0; 0];

    elseif r_k_1 == enum_dyn_model.ct_r || r_k_1 == enum_dyn_model.ct_l
        b = [0; 0; a_cte.x; 0; 0; a_cte.y];
    end

elseif r_k == enum_dyn_model.ct_r || r_k == enum_dyn_model.ct_l
    % Coordinated turn model:
    if r_k_1 == enum_dyn_model.cv
        if r_k == enum_dyn_model.ct_r
            b = [0; 0; 0; 0; -w_cte];
        else
            b = [0; 0; 0; 0;  w_cte];
        end

    elseif r_k_1 == enum_dyn_model.ca
        if r_k == enum_dyn_model.ct_r
            b = [0; 0; 0; 0; -w_cte];
        else
            b = [0; 0; 0; 0;  w_cte];
        end

    elseif r_k_1 == enum_dyn_model.ct_r || r_k_1 == enum_dyn_model.ct_l
        if r_k == enum_dyn_model.ct_r && r_k_1 == enum_dyn_model.ct_l
            b = [0; 0; 0; 0; -w_cte];
        elseif r_k == enum_dyn_model.ct_l && r_k_1 == enum_dyn_model.ct_r
            b = [0; 0; 0; 0;  w_cte];
        else
            b = [0; 0; 0; 0; 0];
        end
    end
end
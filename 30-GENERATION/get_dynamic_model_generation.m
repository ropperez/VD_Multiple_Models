function [F] = get_dynamic_model_generation(r_k,t,w_k,sim_r,ENUM)
% State transition matrix
enum_dyn_model = ENUM.dynamic_model;

switch r_k
    case enum_dyn_model.cv
        % Nearly constant velocity model: 
        if any(sim_r == enum_dyn_model.ct_r) || any(sim_r ==  enum_dyn_model.ct_l)
            F = [1    t    0    0    0; 
                 0    1    0    0    0;
                 0    0    1    t    0;
                 0    0    0    1    0;
                 0    0    0    0    0]; 
        end
    
    case {enum_dyn_model.ct_r, enum_dyn_model.ct_l}
        % Coordinated turn model:  
        % Turn rates
        if r_k == enum_dyn_model.ct_r
            w = -w_k; 
        elseif r_k == enum_dyn_model.ct_l
            w = w_k;
        end

        % Avoiding division by practically zero, w --> 0
        if abs(w) > sqrt(eps)
            % Dynamic model matrix 
            F = [1    sin(w*t)/w        0   -(1-cos(w*t))/w  0; 
                 0    cos(w*t)          0   -sin(w*t)        0;
                 0    (1-cos(w*t))/w    1   sin(w*t)/w       0;
                 0    sin(w*t)          0   cos(w*t)         0;
                 0    0                 0   0                1];    
        else
            % Avoiding division by practically zero, take the limit of w --> 0
            % L'Hopital's rule for the limit
            % Dynamic model matrix 
            F = [1    t    0    0    0; 
                 0    1    0    0    0;
                 0    0    1    t    0;
                 0    0    0    1    0;
                 0    0    0    0    1]; 
        end
end
end

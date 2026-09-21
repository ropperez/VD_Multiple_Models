function [F,b] = get_dynamic_model_static(r_k,t,x_k_1,dim_r_k,sim_r,FILTER_PARAMS,ENUM)

% State transition matrix
enum_dyn_model = ENUM.dynamic_model;

% Initialisation of the output attributes
F = [];
b = 0;

if isempty(x_k_1)
    if r_k == enum_dyn_model.cv
        if any(sim_r == enum_dyn_model.ct_r) || any(sim_r ==  enum_dyn_model.ct_l)
            F = [1    t    0    0    0; 
                 0    1    0    0    0;
                 0    0    1    t    0;
                 0    0    0    1    0;
                 0    0    0    0    0]; 
        end       
    end
else
    if r_k == enum_dyn_model.cv 
        F = FILTER_PARAMS.F_static{r_k};

    elseif r_k == enum_dyn_model.ct_r || r_k == enum_dyn_model.ct_l
        % Coordinated turn model:  
        % Turn rates
        w = x_k_1(dim_r_k);

        if any(sim_r == enum_dyn_model.ca)
            % Avoiding division by practically zero, w --> 0
            if abs(w) > sqrt(eps)
                % Dynamic model matrix 
                sin_w_t = sin(w*t);
                cos_w_t = cos(w*t);
                F = [1    sin_w_t/w        0    0   -(1-cos_w_t)/w   0   0; 
                     0    cos_w_t          0    0   -sin_w_t         0   0;
                     0    0                0    0   0                0   0;
                     0    (1-cos_w_t)/w    0    1   sin_w_t/w        0   0;
                     0    sin_w_t          0    0   cos_w_t          0   0;
                     0    0                0    0   0                0   0;
                     0    0                0    0   0                0   1];    
            else
                % Avoiding division by practically zero, take the limit of w --> 0
                % L'Hopital's rule for the limit
                % Dynamic model matrix 
                F = [1    t    0    0    0    0    0; 
                     0    1    0    0    0    0    0;
                     0    0    0    0    0    0    0;
                     0    0    0    1    t    0    0;
                     0    0    0    0    1    0    0;
                     0    0    0    0    0    0    0;
                     0    0    0    0    0    0    1]; 
            end            
        else
            % Avoiding division by practically zero, w --> 0
            if abs(w) > sqrt(eps)
                % Dynamic model matrix 
                sin_w_t = sin(w*t);
                cos_w_t = cos(w*t);
                F = [1    sin_w_t/w        0   -(1-cos_w_t)/w  0; 
                     0    cos_w_t          0   -sin_w_t        0;
                     0    (1-cos_w_t)/w    1   sin_w_t/w       0;
                     0    sin_w_t          0   cos_w_t         0;
                     0    0                0   0               1];    
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
end
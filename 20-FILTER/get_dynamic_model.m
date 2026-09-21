function [F,b] = get_dynamic_model(r_k,r_k_1,t,x_k,w_cte,FILTER_PARAMS,ENUM)
% F_x: transition matrix for the state
% F_P: transition matrix for the covariance
% b: control parameter

% State transition matrix
enum_dyn_model = ENUM.dynamic_model;

% Initialisation of the output attributes
F = [];
b = 0;

if isempty(r_k_1)
    if r_k == enum_dyn_model.cv 
        F = FILTER_PARAMS.F{r_k,r_k};

    elseif r_k == enum_dyn_model.ct_r || r_k == enum_dyn_model.ct_l
        w_k = x_k(end);

        % Avoiding division by practically zero, w --> 0
        if abs(w_k) > sqrt(eps)
            % Dynamic model matrix
            sin_w_k_t = sin(w_k*t);
            cos_w_k_t = cos(w_k*t);
            F = [1    sin_w_k_t/w_k      0   -(1-cos_w_k_t)/w_k  0 ;
                 0    cos_w_k_t          0   -sin_w_k_t          0 ;
                 0    (1-cos_w_k_t)/w_k  1   sin_w_k_t/w_k       0 ;
                 0    sin_w_k_t          0   cos_w_k_t           0 ;
                 0    0                   0   0                    1];
        else
            % Avoiding division by practically zero, take the limit of w --> 0
            % L'Hopital's rule for the limit
            % Dynamic model matrix
            F = [1    t    0    0    0 ;
                 0    1    0    0    0 ;
                 0    0    1    t    0 ;
                 0    0    0    1    0 ;
                 0    0    0    0    1];
        end
    end

elseif isempty(x_k)
    % Mode in r_k_1
    if r_k == enum_dyn_model.cv
        % Mode in r_k
        F = [1    t    0    0 ;
             0    1    0    0 ;
             0    0    1    t ;
             0    0    0    1];
        
    elseif r_k == enum_dyn_model.ct_r || r_k == enum_dyn_model.ct_l
        % Coordinated turn model:
        F = [1    t    0    0 ;
             0    1    0    0 ;
             0    0    1    t ;
             0    0    0    1 ;
             0    0    0    0];
    end

else

    % Mode in r_k_1
    if r_k == enum_dyn_model.cv
        % Mode in r_k
        if r_k_1 == enum_dyn_model.cv
            F = FILTER_PARAMS.F{r_k,r_k_1};
            b = [0; 0; 0; 0];

        elseif r_k_1 == enum_dyn_model.ct_r || r_k_1 == enum_dyn_model.ct_l
            w_k = x_k(end);

            % Avoiding division by practically zero, w --> 0
            if abs(w_k) > sqrt(eps)
                % Dynamic model matrix
                sin_w_k_t = sin(w_k*t);
                cos_w_k_t = cos(w_k*t);
                F = [1    sin_w_k_t/w_k      0   -(1-cos_w_k_t)/w_k  0 ;
                     0    cos_w_k_t          0   -sin_w_k_t          0 ;
                     0    (1-cos_w_k_t)/w_k  1   sin_w_k_t/w_k       0 ;
                     0    sin_w_k_t          0   cos_w_k_t           0];
            else
                % Avoiding division by practically zero, take the limit of w --> 0
                % L'Hopital's rule for the limit
                % Dynamic model matrix
                F = [1    t    0    0    0 ;
                     0    1    0    0    0 ;
                     0    0    1    t    0 ;
                     0    0    0    1    0];
            end

            b = [0; 0; 0; 0];
        end

    elseif r_k == enum_dyn_model.ct_r || r_k == enum_dyn_model.ct_l
        % Coordinated turn model:
        if r_k_1 == enum_dyn_model.cv
            F = FILTER_PARAMS.F{r_k,r_k_1};

            if r_k == enum_dyn_model.ct_r
                b = [0; 0; 0; 0; -w_cte];
            else
                b = [0; 0; 0; 0;  w_cte];
            end

        elseif r_k_1 == enum_dyn_model.ct_r || r_k_1 == enum_dyn_model.ct_l
            w_k = x_k(end);

            % Avoiding division by practically zero, w --> 0
            if abs(w_k) > sqrt(eps)
                % Dynamic model matrix
                sin_w_k_t = sin(w_k*t);
                cos_w_k_t = cos(w_k*t);
                F = [1    sin_w_k_t/w_k      0   -(1-cos_w_k_t)/w_k  0 ;
                     0    cos_w_k_t          0   -sin_w_k_t          0 ;
                     0    (1-cos_w_k_t)/w_k  1   sin_w_k_t/w_k       0 ;
                     0    sin_w_k_t          0   cos_w_k_t           0 ;
                     0    0                  0   0                   1];
            else
                % Avoiding division by practically zero, take the limit of w --> 0
                % L'Hopital's rule for the limit
                % Dynamic model matrix
                F = [1    t    0    0    0 ;
                     0    1    0    0    0 ;
                     0    0    1    t    0 ;
                     0    0    0    1    0 ;
                     0    0    0    0    1];
            end

            if r_k == enum_dyn_model.ct_r && r_k_1 == enum_dyn_model.ct_l
                b       = [0; 0; 0; 0; -w_cte];
                F(5,5)  = 0;
            elseif r_k == enum_dyn_model.ct_l && r_k_1 == enum_dyn_model.ct_r
                b       = [0; 0; 0; 0;  w_cte];
                F(5,5)  = 0;
            else
                b = [0; 0; 0; 0; 0];
            end
        end
    end
end
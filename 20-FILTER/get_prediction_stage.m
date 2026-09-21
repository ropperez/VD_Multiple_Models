function [x_k_hat, P_k_hat] = get_prediction_stage(x_k_1,P_k_1,r_k,r_k_1,t,dim_r_k_1,FILTER_PARAMS,ENUM,CONSTANTS,SIMULATION)
% Enum and constants
enum_dyn_model  = ENUM.dynamic_model;
dim             = CONSTANTS.dim;
dim_c           = CONSTANTS.dim_c;
w_cte           = SIMULATION.w_gen;

% Posterior in k-1 with correct dimension
x_k_1 = x_k_1(1:dim_r_k_1);
P_k_1 = P_k_1(1:dim_r_k_1,1:dim_r_k_1);  

% Dynamic model transition matrix
[F,b] = get_dynamic_model(r_k,r_k_1,t,x_k_1,w_cte,FILTER_PARAMS,ENUM);

%  Process Noise Matrix
Q = FILTER_PARAMS.Q{r_k,r_k_1};

% State Prediction
x_k_hat = F * x_k_1 + b;

% Jacobian for covariance transformation if constant turn-rate model
if r_k_1 == enum_dyn_model.ct_r || r_k_1 == enum_dyn_model.ct_l
    F_aux = get_jacobian_CT(x_k_1,t,dim,dim_c);
    [m,n] = size(F);
    F     = F_aux(1:m,1:n);  
end

% State Prediction Covariance Matrix
P_k_hat = F * P_k_1 * F' + Q;
end
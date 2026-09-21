function [FILTER] = get_init_filter(x_ini, x_dot_ini, y_ini, y_dot_ini, w_init, sigmaM, dim, ENUM)

% Enums
enum_dyn_model  = ENUM.dynamic_model;

% Initialisation of constants
P_ini_aux       = diag([1 0.1]);
chol_aux        = chol(P_ini_aux)';

% Initial state at k = 1.
FILTER.mean_ini    = [x_ini; x_dot_ini; y_ini; y_dot_ini; w_init]; % state vector definition = [x; x_dot; y; y_dot; w]
chol_ini           = zeros(dim);
chol_ini(1:2,1:2)  = chol_aux;
chol_ini(4:5,4:5)  = chol_aux;

% Initial covariance
FILTER.P_ini       = diag([1 0.1 1 0.1 0.1]);

% Initial dynamic mode selection
if w_init == 0
    FILTER.r_ini       = enum_dyn_model.cv;
elseif w_init < 0
    FILTER.r_ini       = enum_dyn_model.ct_r;
else
    FILTER.r_ini       = enum_dyn_model.ct_l;
end
FILTER.chol_ini    = chol_ini;

% We measure position
FILTER.H_gen       = [1   0   0   0   0 ;
                      0   0   1   0   0];
FILTER.R_k         = (sigmaM)^2*diag([1 1]);

% Generating correlated random numbers
FILTER.chol_R  = chol(FILTER.R_k)';
end
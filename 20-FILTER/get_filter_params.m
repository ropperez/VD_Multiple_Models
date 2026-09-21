function [FILTER_PARAMS] = get_filter_params(sim_r,n_Q,n_t,ENUM,CONSTANTS,SIMULATION)

% Get Constants
dim_nc          = CONSTANTS.dim_nc;
sigmaU          = CONSTANTS.sigmaU;
delta_t         = SIMULATION.delta_t(n_t);

% Init FILTER_PARAMS
FILTER_PARAMS.Q = [];
FILTER_PARAMS.H = [];

% Init attirbutes
F           = cell(max(sim_r),max(sim_r));
F_static    = cell(1,max(sim_r));
Q           = cell(max(sim_r),max(sim_r));
Q_static    = cell(1,max(sim_r));
H           = cell(1,max(sim_r));

% Fill Q cell matrix
for n1 = 1:numel(sim_r)
    r_k             = sim_r(n1);
    H{r_k}          = get_measurement_model(r_k);

    % Dynamic model transition static matrix
    F_static{r_k}   = get_dynamic_model_static(r_k,delta_t,[],[],sim_r,[],ENUM);
    
    %  Process Noise Static Matrix
    Q_static{r_k}   = get_process_noise_static(r_k,sigmaU,delta_t,dim_nc,n_Q,ENUM);

    for n2 = 1:numel(sim_r)
        r_k_1 = sim_r(n2);
        % Dynamic model transition matrix
        F{r_k,r_k_1}  = get_dynamic_model(r_k,r_k_1,delta_t,[],[],[],ENUM);    

        %  Process Noise Matrix
        Q{r_k,r_k_1} = get_process_noise(r_k,r_k_1,sigmaU,delta_t,dim_nc,n_Q,ENUM);        

    end
end

FILTER_PARAMS.F         = F;        % Process noise matrix
FILTER_PARAMS.F_static  = F_static; % Process noise matrix
FILTER_PARAMS.Q         = Q;        % Process noise matrix
FILTER_PARAMS.Q_static  = Q_static; % Process noise matrix
FILTER_PARAMS.H         = H;        % Measurement model matrix

end
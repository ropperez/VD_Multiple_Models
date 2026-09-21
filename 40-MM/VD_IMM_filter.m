function [x_k_mean,P_k_mean,r_k_mean,P_k_mean_all,sum_error2_squared_t] = VD_IMM_filter(n_iter,t,mean_ini,P_ini,x_k_mean,P_k_mean,r_k_mean,P_k_mean_all,X_truth,z_series,R_k,N_r,N_steps,N_iter,sum_error2_squared_t,sim_r,FILTER_PARAMS,ENUM,CONSTANTS,SIMULATION,R_truth,n_Q)
% IMM filter with variable state dimension

% Initialisation
% Constants
dim     = CONSTANTS.dim;
dim_c   = CONSTANTS.dim_c;
dim_z   = CONSTANTS.dim_z;
mu_fil  = CONSTANTS.mu_fil;

% Performance preallocation
x_u_series       = zeros(dim,N_steps,N_iter);
P_u_series       = zeros(dim,dim,N_steps,N_iter);
P_u_trace_series = zeros(dim,dim,N_r,N_steps,N_iter);
r_u_series       = zeros(N_r,N_steps,N_iter);

% State and covariance
x_k     = mean_ini; % r_k = cv -> x = [x; vx; y; vy]
P_k     = P_ini;

% Initialise N_r filters
mu_k_i  = (1/N_r)*ones(N_r,1); % The prior probability for the mode, at the beginning equally probable

% Pre-allocated space
x_k_i           = cell(1,N_r);
P_k_i           = cell(1,N_r);
trace_history   = zeros(N_steps,N_r);

% Initial state
for i = 1:N_r
    r_k                 = sim_r(i);
    [x_k_i{i},P_k_i{i}] = get_init_state_cov(x_k,P_k,r_k,ENUM);  
end

% Save state for performance
x_u_series(:  ,1,n_iter) = x_k;
P_u_series(:,:,1,n_iter) = P_k;
r_u_series(:  ,1,n_iter) = mu_k_i;

for k = 2:N_steps
    % Preallocation of memory for attributes
    z_k             = z_series(:,k);
    x_k_hat_ji      = cell(N_r);
    P_k_hat_ji      = cell(N_r);
    mu_k_hat_i      = zeros(N_r,1);
    mu_k_hat_ji     = zeros(N_r);
    alpha_k_ji      = zeros(N_r);
    x_k_hat         = cell(1,N_r);
    P_k_hat         = cell(1,N_r);
    z_k_hat_i       = zeros(dim_z,N_r);
    z_k_inn_i       = zeros(dim_z,N_r);
    S_k_i           = zeros(dim_z,dim_z,N_r);

    % Save state, coavariance and the prior probability for the mode at k-1 time
    x_k_1_j         = x_k_i;
    P_k_1_j         = P_k_i;
    mu_k_1_j        = mu_k_i;

    % Initialisation of r_k dimension
    dim_r_k_1       = 0;

    % Prediction
    % For each model computes:
    for i = 1:N_r
        r_k = sim_r(i);
        % Dimension of r_k
        dim_r_k     = get_dim_r(r_k);
        max_dim_r   = max(dim_r_k,dim_r_k_1);

        for j = 1:N_r
            r_k_1     = sim_r(j);
            dim_r_k_1 = get_dim_r(r_k_1);

            % Prediction of x_k given r_k, r_k_1
            [x_k_hat_ji{j,i} , P_k_hat_ji{j,i}] = get_prediction_stage(x_k_1_j{j},P_k_1_j{j},r_k,r_k_1,t,dim_r_k_1,FILTER_PARAMS,ENUM,CONSTANTS,SIMULATION);

            % Predicted mode probability
            mu_k_hat_ji(j,i) = mu_fil(j,i)*mu_k_1_j(j);
        end

        % Mixing probability
        mu_k_hat_i(i)   = sum(mu_k_hat_ji(:,i));
        alpha_k_ji(:,i) = mu_k_hat_ji(:,i)/mu_k_hat_i(i);

        % Collapsing the prediction into a Gaussian
        x_k_hat{i} = zeros(dim_r_k,1);
        P_k_hat{i} = zeros(dim_r_k);
        for j = 1:N_r
            x_k_hat{i} = x_k_hat{i} + (alpha_k_ji(j,i) * x_k_hat_ji{j,i});
        end

        for j = 1:N_r
            P_k_hat{i} = P_k_hat{i} + alpha_k_ji(j,i) * (P_k_hat_ji{j,i} + ( (x_k_hat_ji{j,i}-x_k_hat{i}) * (x_k_hat_ji{j,i}-x_k_hat{i})' ) );
        end

        % Prediction measurement
        H_k = FILTER_PARAMS.H{r_k};
        [z_k_hat_i(:,i),S_k_i(:,:,i)] = get_prediction_measurement(x_k_hat{i},P_k_hat{i},H_k,R_k);

        z_k_inn_i(:,i) = z_k-z_k_hat_i(:,i);

        % Update stage
        [x_k_i{i}, P_k_i{i}] = get_update_stage(z_k_inn_i(:,i),S_k_i(:,:,i),x_k_hat{i},P_k_hat{i},H_k,R_k);

        % Trace of Covariance
        trace_history(k,i) = trace(P_k_i{i});
        
        % Save state for performance
        P_u_trace_series(1:dim_r_k,1:dim_r_k,i,k,n_iter) = P_k_i{i};
    end

    % Posterior of the mode
    mu_k_i = div_sum_exp(z_k_inn_i,S_k_i,mu_k_hat_i,N_r);

    % Estimation of the mode
    [~,i_r_k_aux]                    = max(mu_k_i);
    r_k_aux                          = sim_r(i_r_k_aux);
    dim_r_k_aux                      = get_dim_r(r_k_aux);
    x_k                              = zeros(max_dim_r,1);
    P_k                              = zeros(max_dim_r,max_dim_r);
    x_k(1:dim_r_k_aux)               = x_k_i{i_r_k_aux};
    P_k(1:dim_r_k_aux,1:dim_r_k_aux) = P_k_i{i_r_k_aux};

    % Save state for performance
    x_u_series(:  ,k,n_iter) = set_x_k(x_k,r_k_aux,dim);
    P_u_series(:,:,k,n_iter) = P_k;
    r_u_series(:  ,k,n_iter) = mu_k_i;

    % We sum all errors
    sum_error2_squared_t(k)     = sum_error2_squared_t(k)+(X_truth(1,k)-x_k(1))^2+(X_truth(1+dim_c,k)-x_k(1+dim_c))^2;
end
x_k_mean        = x_k_mean     + x_u_series      (:,:    ,n_iter)/N_iter;
P_k_mean        = P_k_mean     + P_u_series      (:,:,:  ,n_iter)/N_iter;
r_k_mean        = r_k_mean     + r_u_series      (:,:    ,n_iter)/N_iter;
P_k_mean_all    = P_k_mean_all + P_u_trace_series(:,:,:,:,n_iter)/N_iter;

end

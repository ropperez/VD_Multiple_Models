function [x_k_mean,P_k_mean,r_k_mean,P_k_mean_all,sum_error2_squared_t] = VD_GPB2_filter(n_iter,t,mean_ini,P_ini,x_k_mean,P_k_mean,r_k_mean,P_k_mean_all,X_truth,z_series,R_k,N_r,N_steps,N_iter,sum_error2_squared_t,sim_r,FILTER_PARAMS,ENUM,CONSTANTS,SIMULATION,R_truth,n_Q)
% GPB2 filter with variable state dimension

% Initialisation
% Constants
mu_fil  = CONSTANTS.mu_fil;
dim     = CONSTANTS.dim;
dim_c   = CONSTANTS.dim_c;
dim_z   = CONSTANTS.dim_z;

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
    z_k                 = z_series(:,k);
    x_k_hat_ji          = cell(N_r);
    P_k_hat_ji          = cell(N_r);
    mu_k_hat_ji         = zeros(N_r);
    x_k_ji              = cell(N_r);
    P_k_ji              = cell(N_r);
    z_k_hat_ji          = cell(N_r);
    z_k_inn_ji          = cell(N_r);
    z_k_inn_j           = zeros(dim_z,N_r);
    S_k_ji              = cell(N_r);
    S_k_j               = zeros(dim_z,dim_z,N_r);

    % Save state, coavariance and the prior probability for the mode at k-1 time
    x_k_1_j         = x_k_i;
    P_k_1_j         = P_k_i;
    mu_k_1_j        = mu_k_i;

    % Initialisation of r_k dimension
    dim_r_k_1       = 0;

    % Prediction
    % Motion model selection
    % r_k = i;

    % For each model computes:
    for i = 1:N_r
        r_k = sim_r(i);
        % Dimension of r_k
        dim_r_k     = get_dim_r(r_k);
        max_dim_r   = max(dim_r_k,dim_r_k_1);

        % r_k_1 = j;
        for j = 1:N_r
            r_k_1       = sim_r(j);
            dim_r_k_1   = get_dim_r(r_k_1);

            % Prediction of x_k given r_k, r_k_1
            [x_k_hat_ji{j,i}, P_k_hat_ji{j,i}] = get_prediction_stage(x_k_1_j{j},P_k_1_j{j},r_k,r_k_1,t,dim_r_k_1,FILTER_PARAMS,ENUM,CONSTANTS,SIMULATION);

            % Prediction measurement
            H_k_i                           = get_measurement_model(r_k);
            [z_k_hat_ji{j,i},S_k_ji{j,i}]   = get_prediction_measurement(x_k_hat_ji{j,i},P_k_hat_ji{j,i},H_k_i,R_k);
            S_k_j(:,:,j)                    = S_k_ji{j,i};

            % Innovation stage
            z_k_inn_ji{j,i} = z_k - z_k_hat_ji{j,i};
            z_k_inn_j(:,j)  = z_k_inn_ji{j,i};

            % Update
            [x_k_ji{j,i}, P_k_ji{j,i}] = get_update_stage(z_k_inn_j(:,j),S_k_j(:,:,j),x_k_hat_ji{j,i},P_k_hat_ji{j,i},H_k_i,R_k);

            % Predicted mode probability
            mu_k_hat_ji(j,i) = mu_fil(j,i)*mu_k_1_j(j);
        end

        % Mixing probability
        w_i = div_sum_exp(z_k_inn_j,S_k_j,mu_k_hat_ji(:,i),N_r);

        % Collapsing step (Gaussian mixture -> in a Gaussian)
        x_k_i{i}    = zeros(dim_r_k,1);
        P_k_i{i}    = zeros(dim_r_k);
    
        % KL divergence minimisation
        x_k_ji_aux  = [x_k_ji{:,i}];
            
        % Mean
        x_k_i{i}    = sum(w_i.' .* x_k_ji_aux,2);
            
        % Covariance
        P_k_ji_aux_1 = zeros(dim_r_k,dim_r_k,N_r);
        for j = 1:N_r
            P_k_ji_aux_1(:,:,j) = (x_k_ji{j,i} - x_k_i{i}) * (x_k_ji{j,i} - x_k_i{i})';
        end
    
        P_k_ji_aux_2 = [P_k_ji{:,i}];
        P_k_ji_aux_2 = reshape(P_k_ji_aux_2,dim_r_k,dim_r_k,N_r);
        
        for j = 1:N_r
            P_k_i{i} = P_k_i{i} + w_i(j) * (P_k_ji_aux_1(:,:,j) + P_k_ji_aux_2(:,:,j));
        end

        trace_history(k,i) = trace(P_k_i{i});
    
        % Save state for performance
        P_u_trace_series(1:dim_r_k,1:dim_r_k,i,k,n_iter) = P_k_i{i};
    end

    % Posterior of the mode
    % Weighted mode
    mu_k_i = sum_div_sumsum_exp(z_k_inn_ji,S_k_ji,mu_k_hat_ji,N_r).';

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

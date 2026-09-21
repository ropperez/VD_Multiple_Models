function [X_truth, R_truth, r_k_truth] = get_ground_truth(t,mean_ini,chol_ini,r_ini,mu_gen,w_gen,sigmaU,sim_r,N_r,dim,dim_c,dim_nc,N_steps,N_q,ENUM)
% Function that generates a ground truth according to the model

% Enums
enum_dyn_model  = ENUM.dynamic_model;

% Preallocation
X_truth             = zeros(dim,N_steps); % True trajectory
R_truth             = zeros(1,N_steps  ); % True mode
r_k_truth           = zeros(N_r,N_steps); % True mode probability

% State at k = 0 
X_truth(:,1)        = mean_ini(:)+chol_ini*randn(dim,1);
R_truth(:,1)        = r_ini;
r_k_truth(r_ini,1)  = 1;

% Dynamic model at k = 0
r_k = r_ini;

% Target and mode truth
for k = 2:N_steps
    % Mode at k-1
    r_k_1       = r_k;
    r_k_1_aux   = sim_r(r_k_1);

    % Motion model selection
    r_k     = get_mode(r_k_1,mu_gen);
    r_k_aux = sim_r(r_k);

    %  Dynamic function (matrix) of the state.
    F_k     = get_dynamic_model_generation(r_k_aux,t,w_gen,sim_r,ENUM);

    % Process noise, covariance matrix
    Q       = get_process_noise_generation(r_k_aux,sigmaU,t,dim_nc,sim_r,N_q,ENUM);

    % Choleski decomposition 
    if r_k_aux == enum_dyn_model.cv
        chol_Q          = zeros(dim);
        Q_aux           = Q(1:2,1:2);
        chol_Q_aux      = chol(Q_aux)';
        chol_Q(1:2,1:2) = chol_Q_aux;
        chol_Q(3:4,3:4) = chol_Q_aux;
    elseif r_k_aux == enum_dyn_model.ct_r || r_k_aux == enum_dyn_model.ct_l
        chol_Q = chol(Q)';
    else
        error('Dynamic model not supported for generation.')
    end

    % Transition matrix
    X_truth(:,k) = F_k*X_truth(:,k-1);

    % Process noise
    Q_truth      = chol_Q*randn(dim,1);

    % Dynamic model
    X_truth(:,k) = X_truth(:,k)+Q_truth;
    if r_k_aux == enum_dyn_model.ct_r && r_k_1_aux ~= enum_dyn_model.ct_r
        X_truth(dim,k) = - w_gen + Q_truth(dim);
    elseif r_k_aux == enum_dyn_model.ct_l && r_k_1_aux ~= enum_dyn_model.ct_l
        X_truth(dim,k) =   w_gen + Q_truth(dim);
    end

    % Measurement noise matrix
    R_truth(:,k) = r_k_aux;

    % Mode
    if r_k_aux > N_r
        r_k_truth(N_r,k) = 1;
    else
        r_k_truth(r_k_aux,k) = 1;
    end
end

% Figure 1
% Target truth representation for paper
figure('Name','Target Truth');
% Position x, y
box on, grid on, axis equal, hold on
title('Position')
xlabel('P_{x} [m]'),ylabel('P_{y} [m]')

% Trajectory Mode
legend_tag = ({'CV','CT_R','CT_L', 'CA'});
legend_tag = legend_tag(sim_r);

% Position x, y
for n = 1:numel(sim_r)
    i_r_k                   = R_truth(1,:) == sim_r(n);
    X_truth_aux             = nan(2,numel(R_truth(1,:)));
    X_truth_aux(:,i_r_k)    = X_truth([1 1+dim_c],i_r_k);
    scatter(X_truth_aux(1,:),X_truth_aux(2,:),"filled")
end

legend(legend_tag)
drawnow

end
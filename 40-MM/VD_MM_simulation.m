function VD_MM_simulation(sim_r,sim_filter,tag_filter,ENUM)

% Enums
enum_filter     = ENUM.enum_filter;

% Get simulation params
SIMULATION    = get_simulation_params() ;
N_steps       = SIMULATION.N_steps      ; % Number of steps for each iteration
N_iter        = SIMULATION.N_iter       ; % Number of iterations for each iteration
delta_t       = SIMULATION.delta_t      ; % Sampling time
x_ini         = SIMULATION.x_ini        ; % Initial target x position           [m]
x_dot_ini     = SIMULATION.x_dot_ini    ; % Initial target x velocity           [m/s]
y_ini         = SIMULATION.y_ini        ; % Initial target y position           [m]
y_dot_ini     = SIMULATION.y_dot_ini    ; % Initial target y velocity           [m/s]
w_init        = SIMULATION.w_init       ; % Initial target turn rate            [rad/s]
w_gen         = SIMULATION.w_gen        ; % Constant target turn rate           [rad/s];         Constant for b control vector in motion model [rad/s]

N_r           = numel(sim_r)            ; % Number of dynamic models for update

% Get constants
CONSTANTS   = get_init_constants();
dim         = CONSTANTS.dim      ; % Number of state components
dim_c       = CONSTANTS.dim_c    ; % Number of state cartesian components
dim_nc      = CONSTANTS.dim_nc   ; % Number of cartesian components
sigmaU      = CONSTANTS.sigmaU   ; % Standard deviation of the generated noise
sigmaM      = CONSTANTS.sigmaM   ; % Standard deviation of the measurement noise
mu_gen      = CONSTANTS.mu_gen   ; % Mode transition matrix p(r_k = i/r_{k-1} = j) for generation of target truth
gamma       = CONSTANTS.gamma    ; % Scaling factor for the process noise
dim_t       = numel(delta_t)     ; % Number of sample time simulated
dim_M       = numel(sigmaM)      ; % Number of measurement noise simulated
dim_Q       = numel(gamma)       ; % Number of process noise simulated

% Number of sampling time
for n_delta_t = 1:dim_t
    % Simulation status
    fprintf('Sampling time %d of %d \n',n_delta_t,dim_t);

    % Number of measuremtn noise
    for n_M = 1:dim_M
        % Simulation status
        fprintf('    Measurement Noise %d of %d \n',n_M,dim_M);

        % Filter initialisation
        FILTER      = get_init_filter(x_ini, x_dot_ini, y_ini, y_dot_ini, w_init, sigmaM(n_M), dim, ENUM);
        mean_ini    = FILTER.mean_ini;
        P_ini       = FILTER.P_ini   ;
        r_ini       = FILTER.r_ini   ;
        chol_ini    = FILTER.chol_ini;
        H_gen       = FILTER.H_gen   ;
        R_k         = FILTER.R_k     ;
        chol_R      = FILTER.chol_R  ;
    
        % Number of process noise
        for n_Q = 1:dim_Q
            % Fix seed for simulation of random numbers
            rng('default')
            rng(20)

            % Simulation status
            fprintf('        Process Noise %d of %d \n',n_Q,dim_Q);
            % We generate ground truth according to the model
            [X_truth, R_truth, ~] = get_ground_truth(delta_t(n_delta_t),mean_ini,chol_ini,r_ini,mu_gen,w_gen,sigmaU,sim_r,N_r,dim,dim_c,dim_nc,N_steps,n_Q,ENUM);
            
            % Get measurements for filtering
            z_series = get_measurements_generation(X_truth,chol_R,H_gen,N_steps,N_iter);

            % Get filter params
            FILTER_PARAMS = get_filter_params(sim_r,n_Q,n_delta_t,ENUM,CONSTANTS,SIMULATION);
                      
            % MC Simulation
            for n_sim_filter = 1:numel(sim_filter)
                % Simulation status
                fprintf(['            Filter %d of %d -> ' tag_filter{sim_filter(n_sim_filter)} '\n'],n_sim_filter,numel(sim_filter));
            
                % Position error (squared)
                sum_error2_squared_t = zeros(N_steps,1);
            
                % Performance
                x_k_mean         = 0; % State mean for MC
                P_k_mean         = 0; % Covariance mean for MC
                r_k_mean         = 0; % Mode mean for MC
                P_k_mean_all     = 0; % Covariance mean trace for MC
                            
                % Iterations
                for n_iter = 1:N_iter
                    % MM Filter selection
                    switch sim_filter(n_sim_filter)
                        % VD-IMM filter
                        case enum_filter.vd_imm
                            [x_k_mean,P_k_mean,r_k_mean,P_k_mean_all,sum_error2_squared_t] = VD_IMM_filter (n_iter,delta_t(n_delta_t),mean_ini,P_ini,x_k_mean,P_k_mean,r_k_mean,P_k_mean_all,X_truth,z_series(:,:,n_iter),R_k,N_r,N_steps,N_iter,sum_error2_squared_t,sim_r,FILTER_PARAMS,ENUM,CONSTANTS,SIMULATION,R_truth,n_Q);
                            plot_n = 1;
                        % VD-GPB2 filter
                        case enum_filter.vd_gpb2
                            [x_k_mean,P_k_mean,r_k_mean,P_k_mean_all,sum_error2_squared_t] = VD_GPB2_filter(n_iter,delta_t(n_delta_t),mean_ini,P_ini,x_k_mean,P_k_mean,r_k_mean,P_k_mean_all,X_truth,z_series(:,:,n_iter),R_k,N_r,N_steps,N_iter,sum_error2_squared_t,sim_r,FILTER_PARAMS,ENUM,CONSTANTS,SIMULATION,R_truth,n_Q);
                            plot_n = 2;
                        otherwise
                            error('Filter not supported yet.')
                    end
                end
            
                % Compute trace of filter covariance
                % Number of modes adjustment in static filters
                N_r_aux = size(P_k_mean_all,3);
                P_k_mean_trace = zeros(1,N_steps,N_r_aux);
                for n_steps = 1:N_steps
                    for n_r = 1:N_r_aux
                        P_k_mean_trace(:,n_steps,n_r) = trace(P_k_mean_all(:,:,n_r,n_steps));
                    end
                end
                
                % Save RMSE
                rmse_error_t   = sqrt(sum_error2_squared_t/N_iter);
                
                % RMS and Trajectory mode representation
                % Set color order
                color_order = ['g' 'm'];

                if plot_n == 1
                    f_rms_mode = figure('Name','RMS Position Error - Real Mode');
                    tl_rms_mode = tiledlayout(f_rms_mode,2,1);
    
                    % RMS Position
                    s1_rms_mode = nexttile(tl_rms_mode);
                    box on, grid on, grid minor, hold on
                    title(['$T = $ ' num2str(delta_t(n_delta_t)) ' s; $\sigma_m = $ ' num2str(sigmaM(n_M)) ' m; $\gamma = $ ' num2str(gamma(n_Q))],'Interpreter','latex')
                    xlabel('Time step'), ylabel('RMSE_k (m)')   
                    xlim([1 100])
                    yscale log
        
                    % Trajectory Mode
                    s2_rms_mode = nexttile(tl_rms_mode);
                    box on, grid on, grid minor, hold on            
                    title('Trajectory Mode','Interpreter','latex')
                    xlabel('Time Step')
                    yticks([1 2 3 4])
                    yticklabels({'CV','CT\_R','CT\_L', 'CA'})
                    ylim([min(sim_r)-0.5 max(sim_r)+0.5])
                end
                    % Trajectory Mode
                scatter(s2_rms_mode,1:N_steps,squeeze(R_truth),'b.')
        
                % RMS representation
                if plot_n == 1
                    f_rms = figure('Name','RMS Position Error');
                    f_rms_ax      = axes('Parent',f_rms);
                    yscale log
    
                    % RMS Position
                    box on, grid on, grid minor, hold on
                    title(['$T = $ ' num2str(delta_t(n_delta_t)) ' s; $\sigma_m = $ ' num2str(sigmaM(n_M)) ' m; $\gamma = $ ' num2str(gamma(n_Q))],'Interpreter','latex')
                    xlabel('Time step'), ylabel('RMSE_k (m)')                          
                end

                % Real Mode representation
                if plot_n == 1
                    f_mode = figure('Name','Real Mode');
                    f_mode_ax      = axes('Parent',f_mode);

                    % Trajectory Mode
                    box on, grid on, grid minor, hold on
                    title('Trajectory Mode')
                    xlabel('Time Step')
                    yticks([1 2 3 4])
                    yticklabels({'CV','CT\_R','CT\_L'})
                    ylim([min(sim_r)-0.2 max(sim_r)+0.2])
                end

                % Trajectory Mode
                scatter(f_mode_ax,1:N_steps,squeeze(R_truth),'b.')     

                plot(s1_rms_mode,rmse_error_t,[color_order(sim_filter(n_sim_filter)) '.-'])
                plot(f_rms_ax   ,rmse_error_t,[color_order(sim_filter(n_sim_filter)) '.-'])
            end
                                           
            % Set legend for the filter name
            legend_tag_filter = {'VD-IMM';'VD-GPB2'};

            legend(s1_rms_mode,legend_tag_filter(sim_filter),'location','best');
            legend(f_rms_ax   ,legend_tag_filter(sim_filter),'location','best');
        end
    end
end
end
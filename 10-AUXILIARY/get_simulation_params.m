function [SIMULATION] = get_simulation_params()
% Initialisation of parameters

% MC parameters
SIMULATION.N_steps      = 100       ; % Number of steps for each iteration
SIMULATION.N_iter       = 100       ; % Number of iterations for each trajectory
SIMULATION.delta_t      = 0.5       ; % Sampling time, it could be a vector [0.5 1] [s]

% Simulation generation turn rate
SIMULATION.w_gen        = 15 *pi/180; % turn rate [rad/s]

% Trajectory initialisation
SIMULATION.x_ini        = 0         ; % Initial target x position [m]     
SIMULATION.x_dot_ini    = 50        ; % Initial target x velocity [m/s]    
SIMULATION.y_ini        = 0         ; % Initial target y position [m]    
SIMULATION.y_dot_ini    = 50        ; % Initial target y velocity [m/s]    
SIMULATION.w_init       = 0         ; % Initial target turn rate  [rad/s]
end
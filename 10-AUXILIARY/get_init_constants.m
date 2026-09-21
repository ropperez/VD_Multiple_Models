function [CONSTANTS] = get_init_constants()

% Initialisation of constants
CONSTANTS.dim_c  = 2; % Number of state cartesian components [x; x_dot]
CONSTANTS.dim_nc = 2; % Number of cartesian components       [x; y]; 
CONSTANTS.dim_z  = 2; % Number of measurements components    [x; y]; 
CONSTANTS.dim    = 5; % Number of state components           [x; x_dot; y; y_dot; w]

% Scaling factor for process noise
CONSTANTS.gamma         = 1e-2; % it could be a vector, [1e-2     10^-1.75    10^-1.5     10^-1.25    10^-1   10^-0.75    10^-0.5     10^-0.25    1    5];

% Process Noise
CONSTANTS.sigmaU.cv     = CONSTANTS.gamma * 1;         % Standard deviation of the generated noise        [m/s^3/2]
CONSTANTS.sigmaU.ct_r   = CONSTANTS.gamma * 1*pi/180;  % Standard deviation of the generated noise, right [rad/s^1/2]
CONSTANTS.sigmaU.ct_l   = CONSTANTS.gamma * 1*pi/180;  % Standard deviation of the generated noise, left  [rad/s^1/2]

% Measurement noise, it could be a vector
CONSTANTS.sigmaM        = 0.5;                         % Standard deviation of the measurement noise, it could be a vector [0.5  1] [m]

% Mode transition probability matrix for generation    
CONSTANTS.mu_gen        = [0.9      0.05    0.05 ;     % Mode transition matrix p(r_k = i/r_{k-1} = j) for generation
                           0.05     0.9     0.05 ;
                           0.05     0.05    0.9 ];

% Mode transition probability matrix for filtering  
CONSTANTS.mu_fil        = [0.9      0.05    0.05 ;     % Mode transition matrix p(r_k = i/r_{k-1} = j) for filtering
                           0.05     0.9     0.05 ;
                           0.05     0.05    0.9 ];  
end
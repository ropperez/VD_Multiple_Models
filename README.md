# Variable Dimension IMM and GPB2 Filters for Tracking in Multiple Model Systems

This repository contains the official MATLAB implementation for the paper:  
**"Variable Dimension IMM and GPB2 Filters for Tracking in Multiple Model Systems"**  
Published in **IEEE Transactions on Aerospace and Electronic Systems**.
* **Paper Authors:** Roberto Pérez-Pérez and Ángel F. García-Fernández
* **DOI:** https://doi.org/10.1109/TAES.2026.3735979
* **Code Implementation & Maintainer:** Roberto Pérez-Pérez  
* **Affiliation:** Information Processing and Telecommunications Center, ETSI de Telecomunicación, Universidad Politécnica de Madrid, 28040 Madrid, Spain.  
* **Repository URL:** [https://github.com/ropperez/VD_Multiple_Models](https://github.com/ropperez/VD_Multiple_Models)

---

## 📌 Abstract

This paper presents mathematically principled filters for multiple-model systems with states of different dimensionality, specifically the variable-dimension interacting multiple model (**VD-IMM**) filter and the variable-dimension generalised pseudo-Bayesian filter of order 2 (**VD-GPB2**). 

To do so, we first provide a Bayesian modelling of a variable dimensional dynamic system and its measurements. Then, for variable dimensional linear Gaussian dynamic and measurement models:
* The **VD-IMM filter** is derived by assuming a posterior that has a Gaussian density for each mode, and then performing a Kullback-Leibler Divergence (KLD) minimisation after each prediction step to keep the Gaussian density form for each mode. Subsequently, the Bayesian update step is performed, keeping a Gaussian density for each mode.
* The **VD-GPB2 filter** considers the same model as the VD-IMM filter and also assumes that the posterior for each mode is Gaussian. In contrast, the VD-GPB2 filter propagates a Gaussian mixture for each mode in the prediction step. In the update step, the VD-GPB2 filter performs a KLD minimisation to have a Gaussian density for each mode.

Simulation results demonstrate the advantages of the VD-IMM and VD-GPB2 filters over previous state-of-the-art alternatives.

---

## 🛠️ Requirements & Installation

### Environment
* **MATLAB R2025b** (or compatible release)

### Repository Setup
Clone this repository to your local machine:
```bash
git clone https://github.com/ropperez/VD_Multiple_Models.git
cd VD_Multiple_Models
```

---

## 📁 Repository Structure

```text
VD_Multiple_Models/
├── 10-AUXILIARY/           # Enumerations, initial constants, and simulation parameters
├── 20-FILTER/              # Filter components (prediction, update, dynamic and measurement models)
├── 30-GENERATION/          # Trajectory ground truth and measurement generation
├── 40-MM/                  # Multiple-Model filter routines (VD-IMM, VD-GPB2, MM simulation runner)
├── 50-UTILS/               # Mathematical utilities and log-likelihood calculation
└── test_VD_MM_simulation.m # Main test script to configure and launch simulations
```

---

## ⚙️ Simulation Configuration & Running
### Simulation Configuration
You can configure the active filters and the dynamic model modes directly inside `test_VD_MM_simulation.m`:

```matlab
% Filter type simulation selection
sim_filter = [enum_filter.vd_imm enum_filter.vd_gpb2];

% Models to be generated
sim_r = [enum_dyn_model.cv enum_dyn_model.ct_r enum_dyn_model.ct_l];
```
#### Filter and Modes Description
- `sim_filter`: Selects the filtering algorithms to be evaluated and compared in the simulation run:

  - `enum_filter.vd_imm`: VD-IMM filter.

  - `enum_filter.vd_gpb2`: VD-GPB2 filter.

- `sim_r` : Specifies the dynamic model modes used to synthesize the target's trajectory ground truth:

  - `enum_dyn_model.cv`: Constant Velocity model.

  - `enum_dyn_model.ct_r`: Coordinated Turn Right model.

  - `enum_dyn_model.ct_l`: Coordinated Turn Left model.

#### Simulation Parameters Description
Simulation parameters are configured in `10-AUXILIARY/get_simulation_params.m`:

```matlab
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
```

##### Parameter Description
- `SIMULATION.N_steps`: Number of time steps evaluated during each Monte Carlo run.
- `SIMULATION.N_iter`: Total number of independent Monte Carlo trials for performance evaluation.
- `SIMULATION.delta_t`: Sampling time ($T$) in seconds. Note: A vector of values (e.g., [0.5, 1.0]) can be provided to evaluate and generate simulation results across multiple sampling intervals.
- `SIMULATION.w_gen`: Turn rate (in rad/s) applied during maneuvering trajectory phases.
- `SIMULATION.x_ini` / `SIMULATION.y_ini`: Initial target position coordinates ($x, y$) in meters.
- `SIMULATION.x_dot_ini` / `SIMULATION.y_dot_ini`: Initial target velocity components ($\dot{x}, \dot{y}$) in m/s.
- `SIMULATION.w_init`: Initial target turn rate ($\omega$) in rad/s.

#### Constant Parameters Description
Constant parameters are configured in `10-AUXILIARY/get_init_constants.m`:

```matlab
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
```

##### Parameter Description
- `CONSTANTS.gamma`: Scaling factor for the process noise. A vector of values can be passed (e.g., [1e-2, 10^-1.5, 1]), which will trigger and execute an independent Monte Carlo simulation run for each value in the vector.
- `CONSTANTS.sigmaM`: Standard deviation of the measurement noise ($m$). Similar to process noise, a vector with multiple noise levels can be specified (e.g., [0.5, 1]), running a Monte Carlo simulation for each measurement noise parameter.
- `CONSTANTS.mu_gen`: Mode transition probability matrix used for ground-truth trajectory generation ($p(r_k = i \mid r_{k-1} = j)$). Fully configurable to model target maneuver switching frequencies.
- `CONSTANTS.mu_fil`: Mode transition probability matrix utilized internally by the estimation filters ($p(r_k = i \mid r_{k-1} = j)$). Configurable to evaluate filter performance under potential model mismatches relative to `mu_gen`.
##### State Dimensions:
- `CONSTANTS.dim_c`: Cartesian components per dimension ($[x, \dot{x}]$).
- `CONSTANTS.dim_z`: Measured components ($[x, y]$).
- `CONSTANTS.dim`: Total dimension of the extended state ($[x, \dot{x}, y, \dot{y}, \omega]$).

### Running
To run the simulation:
1. Open MATLAB and navigate to the repository root directory (`VD_Multiple_Models`).
2. Run the main simulation test script from the Command Window:
```matlab
test_VD_MM_simulation
```
---
## 📊 Generated Output Figures

Upon executing the simulation test script (`test_VD_MM_simulation.m`), the framework automatically generates and plots three key performance figures:

### 1. RMS Position Error & Real Mode (`RMS Position Error - Real Mode`)
A combined tiled layout (`2x1`) illustrating both estimation accuracy and target maneuver switches over time:
* **Top Subplot (RMSE):** Displays the root-mean-square error ($\text{RMSE}_k$) in position (meters) across time steps on a logarithmic scale for the evaluated filters. The plot title reflects the specific scenario parameters: sampling interval ($T$), measurement noise ($\sigma_m$), and process noise scaling ($\gamma$).
* **Bottom Subplot (Trajectory Mode):** Displays the ground-truth target mode sequence ($r_k$) across time steps, mapped to dynamic behaviors: `CV` (Constant Velocity), `CT_R` (Coordinated Turn Right) and `CT_L` (Coordinated Turn Left).

### 2. RMS Position Error (`RMS Position Error`)
A standalone logarithmic plot dedicated exclusively to comparing the position tracking accuracy ($\text{RMSE}_k$) among the active multiple-model filters (`VD-IMM`, `VD-GPB2`, etc.).

### 3. Real Trajectory Mode (`Real Mode`)
A standalone scatter plot showing the active ground-truth motion model for each time step of the generated trajectory.

---

## 📄 Citation & Paper Link

If you use this code or reference our work in your research, please cite our paper:

```bibtex
@article{perez2026variable,
  title={Variable Dimension IMM and GPB2 Filters for Tracking in Multiple Model Systems},
  author={P{\'e}rez-P{\'e}rez, Roberto and Garc{\'i}a-Fern{\'a}ndez, {\'A}ngel F.},
  journal={IEEE Transactions on Aerospace and Electronic Systems},
  year={2026},
  doi={10.1109/TAES.2026.3735979}
}
```

🔗 **IEEE Xplore / Official DOI:** [https://doi.org/10.1109/TAES.2026.3735979](https://doi.org/10.1109/TAES.2026.3735979)

---


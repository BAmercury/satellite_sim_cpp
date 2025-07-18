function [A, B, C, D, state_names, input_names, system_info] = build_multi_tank_slosh_state_space(spacecraft_params, tank_params)
% BUILD_MULTI_TANK_SLOSH_STATE_SPACE Creates state space model for spacecraft with multiple propellant tanks
%
% Inputs:
%   spacecraft_params: Structure with fields:
%     - I: 3x3 inertia matrix [kg*m^2]
%     - M0: Total spacecraft mass [kg]
%     - g: Effective acceleration (thrust/mass) [m/s^2]
%
%   tank_params: Structure array (one element per tank) with fields:
%     - m: Slosh mass [kg]
%     - l: Pendulum length [m]
%     - xH, yH, zH: Hinge location w.r.t. spacecraft CG [m]
%     - zeta: Viscous damping factor [dimensionless]
%     - name: String identifier for the tank (optional)
%
% Outputs:
%   A, B, C, D: State space matrices
%   state_names: Cell array of state variable names
%   input_names: Cell array of input names
%   system_info: Structure with system information
%
% State vector: [psi phi theta wx wy wz psi_p1 theta_p1 psi_p1_dot theta_p1_dot … psi_pN theta_pN psi_pN_dot theta_pN_dot]
% Input vector: [Tcx Tcy Tcz fx_thr fy_thr fz_thr]

% Extract spacecraft parameters
I_matrix = spacecraft_params.I;  % Full 3x3 inertia matrix
M0 = spacecraft_params.M0;
g = spacecraft_params.g;

% Calculate inverse inertia matrix (CORRECT approach for non-diagonal inertia)
I_inv = inv(I_matrix);
inv_Ix = I_inv(1,1);  % 1/I_effective_x (accounts for cross-coupling)
inv_Iy = I_inv(2,2);  % 1/I_effective_y
inv_Iz = I_inv(3,3);  % 1/I_effective_z

% Number of tanks
n_tanks = length(tank_params);

% Pre-calculate tank properties
tank_info = struct();
for i = 1:n_tanks
% Extract tank parameters
tank_info(i).m = tank_params(i).m;
tank_info(i).l = tank_params(i).l;
tank_info(i).xH = tank_params(i).xH;
tank_info(i).yH = tank_params(i).yH;
tank_info(i).zH = tank_params(i).zH;
tank_info(i).zeta = tank_params(i).zeta;

```
% Calculate derived parameters for each tank
tank_info(i).omega_sl = sqrt(g/tank_info(i).l * (1 + tank_info(i).m/M0));

% Acceleration sensitivities (using correct inverse inertia elements)
tank_info(i).mu_P_accn_x = tank_info(i).m * tank_info(i).yH * g * inv_Ix;
tank_info(i).mu_P_accn_z = tank_info(i).m * tank_info(i).yH * g * inv_Iz;

% Rate sensitivities
tank_info(i).mu_P_rate_x = tank_info(i).mu_P_accn_x / tank_info(i).omega_sl;
tank_info(i).mu_P_rate_z = tank_info(i).mu_P_accn_z / tank_info(i).omega_sl;

% Tank name for identification
if isfield(tank_params(i), 'name')
    tank_info(i).name = tank_params(i).name;
else
    tank_info(i).name = sprintf('Tank%d', i);
end
```

end

% State vector dimension: 6 rigid body + 4*n_tanks pendulum states
N = 6 + 4*n_tanks;

% Initialize matrices
A = zeros(N, N);
B = zeros(N, 6);
C = eye(N);  % Full state output
D = zeros(N, 6);

%% Build A matrix

% Attitude kinematics: [psi; phi; theta] = [wx; wy; wz]
A(1:3, 4:6) = eye(3);

% Slosh coupling to spacecraft dynamics
for i = 1:n_tanks
% State indices for tank i
psi_p_idx = 6 + (i-1)*4 + 1;    % psi_p for tank i
theta_p_idx = 6 + (i-1)*4 + 2;  % theta_p for tank i

```
% Effect of psi_p on wx (yaw rate) - equation from paper
A(4, psi_p_idx) = A(4, psi_p_idx) - tank_info(i).mu_P_accn_x;

% Effect of theta_p on wz (pitch rate) - equation from paper  
A(6, theta_p_idx) = A(6, theta_p_idx) - tank_info(i).mu_P_accn_z;
```

end

% Pendulum dynamics for each tank
for i = 1:n_tanks
% State indices for tank i
base_idx = 6 + (i-1)*4;
psi_p_idx = base_idx + 1;       % psi_p
theta_p_idx = base_idx + 2;     % theta_p  
psi_p_dot_idx = base_idx + 3;   % psi_p_dot
theta_p_dot_idx = base_idx + 4; % theta_p_dot

```
% Pendulum kinematics: [psi_p; theta_p] = [psi_p_dot; theta_p_dot]
A(psi_p_idx, psi_p_dot_idx) = 1;
A(theta_p_idx, theta_p_dot_idx) = 1;

% Pendulum damping
A(psi_p_dot_idx, psi_p_dot_idx) = -2*tank_info(i).zeta*tank_info(i).omega_sl;
A(theta_p_dot_idx, theta_p_dot_idx) = -2*tank_info(i).zeta*tank_info(i).omega_sl;

% Pendulum restoring forces with coupling from all tanks
% For psi_p equation (yaw-roll plane) - diagonal term
A(psi_p_dot_idx, psi_p_idx) = -(1/tank_info(i).l) * ...
    ((1 + tank_info(i).m/M0)*g + (tank_info(i).yH - tank_info(i).l)*tank_info(i).mu_P_accn_z);

% For theta_p equation (roll-pitch plane) - diagonal term
A(theta_p_dot_idx, theta_p_idx) = -(1/tank_info(i).l) * ...
    ((1 + tank_info(i).m/M0)*g + (tank_info(i).yH - tank_info(i).l)*tank_info(i).mu_P_accn_x);

% Cross-coupling terms with other tanks (from equations 45-46)
for j = 1:n_tanks
    if i ~= j
        other_psi_p_idx = 6 + (j-1)*4 + 1;
        other_theta_p_idx = 6 + (j-1)*4 + 2;
        
        % Cross-coupling in psi equation
        A(psi_p_dot_idx, other_psi_p_idx) = -(1/tank_info(i).l) * ...
            (tank_info(j).m/M0 * g + (tank_info(i).yH - tank_info(i).l)*tank_info(j).mu_P_accn_z);
        
        % Cross-coupling in theta equation  
        A(theta_p_dot_idx, other_theta_p_idx) = -(1/tank_info(i).l) * ...
            (tank_info(j).m/M0 * g + (tank_info(i).yH - tank_info(i).l)*tank_info(j).mu_P_accn_x);
    end
end
```

end

%% Build B matrix

% Control torques on spacecraft attitude rates (using inverse inertia elements)
B(4, 1) = inv_Ix;  % Tcx affects wx
B(5, 2) = inv_Iy;  % Tcy affects wy  
B(6, 3) = inv_Iz;  % Tcz affects wz

% Thruster forces affecting pendulum dynamics for each tank
for i = 1:n_tanks
base_idx = 6 + (i-1)*4;
psi_p_dot_idx = base_idx + 3;   % psi_p_dot
theta_p_dot_idx = base_idx + 4; % theta_p_dot

```
% psi_p equation (from equation 47 in paper)
B(psi_p_dot_idx, 2) = -tank_info(i).zH/(tank_info(i).l * (1/inv_Iy));  % fy_thr effect
B(psi_p_dot_idx, 3) = (tank_info(i).yH - tank_info(i).l)/(tank_info(i).l * (1/inv_Iz));  % fz_thr effect
B(psi_p_dot_idx, 4) = -1/(M0*tank_info(i).l);  % fx_thr effect

% theta_p equation  
B(theta_p_dot_idx, 1) = (tank_info(i).yH - tank_info(i).l)/(tank_info(i).l * (1/inv_Ix));  % fx_thr effect
B(theta_p_dot_idx, 2) = -tank_info(i).xH/(tank_info(i).l * (1/inv_Iy));  % fy_thr effect
B(theta_p_dot_idx, 6) = 1/(M0*tank_info(i).l);  % fz_thr effect
```

end

%% Generate state and input names
state_names = {‘psi’, ‘phi’, ‘theta’, ‘wx’, ‘wy’, ‘wz’};

for i = 1:n_tanks
tank_name = tank_info(i).name;
state_names{end+1} = sprintf(‘psi_p_%s’, tank_name);
state_names{end+1} = sprintf(‘theta_p_%s’, tank_name);
state_names{end+1} = sprintf(‘psi_p_dot_%s’, tank_name);
state_names{end+1} = sprintf(‘theta_p_dot_%s’, tank_name);
end

input_names = {‘Tcx’, ‘Tcy’, ‘Tcz’, ‘fx_thr’, ‘fy_thr’, ‘fz_thr’};

%% Package system information
system_info.n_tanks = n_tanks;
system_info.tank_info = tank_info;
system_info.n_states = N;
system_info.inertia_matrix = I_matrix;
system_info.inverse_inertia = I_inv;
system_info.spacecraft_mass = M0;
system_info.total_slosh_mass = sum([tank_info.m]);

% Check if inertia matrix is diagonal
system_info.is_diagonal_inertia = all(all(abs(I_matrix - diag(diag(I_matrix))) < 1e-10));
if ~system_info.is_diagonal_inertia
fprintf(‘Note: Non-diagonal inertia matrix detected. Using proper inverse inertia approach.\n’);
end

end

%% Example usage for multiple tanks
function example_multi_tank_usage()
% Define spacecraft with non-diagonal inertia matrix
spacecraft_params.I = [1000  50   20;   % Non-diagonal inertia matrix
50   2000  30;
20   30   1500];
spacecraft_params.M0 = 2000;    % kg
spacecraft_params.g = 9.81;     % m/s^2

```
% Define 4 tanks: 2 fuel tanks + 2 oxidizer tanks
tank_params(1).m = 100;     % Fuel tank 1
tank_params(1).l = 0.6;
tank_params(1).xH = 0.5;
tank_params(1).yH = 1.2;
tank_params(1).zH = 0.3;
tank_params(1).zeta = 0.03;
tank_params(1).name = 'Fuel1';

tank_params(2).m = 100;     % Fuel tank 2 (symmetric)
tank_params(2).l = 0.6;
tank_params(2).xH = -0.5;
tank_params(2).yH = 1.2;
tank_params(2).zH = 0.3;
tank_params(2).zeta = 0.03;
tank_params(2).name = 'Fuel2';

tank_params(3).m = 150;     % Oxidizer tank 1 (different properties)
tank_params(3).l = 0.8;
tank_params(3).xH = 0.4;
tank_params(3).yH = 2.0;
tank_params(3).zH = -0.2;
tank_params(3).zeta = 0.05;
tank_params(3).name = 'Ox1';

tank_params(4).m = 150;     % Oxidizer tank 2
tank_params(4).l = 0.8;
tank_params(4).xH = -0.4;
tank_params(4).yH = 2.0;
tank_params(4).zH = -0.2;
tank_params(4).zeta = 0.05;
tank_params(4).name = 'Ox2';

% Build multi-tank state space model
[A, B, C, D, state_names, input_names, system_info] = ...
    build_multi_tank_slosh_state_space(spacecraft_params, tank_params);

% Display results
fprintf('Multi-Tank Slosh System Built Successfully!\n');
fprintf('Number of tanks: %d\n', system_info.n_tanks);
fprintf('Total states: %d\n', system_info.n_states);
fprintf('Total slosh mass: %.1f kg (%.1f%% of spacecraft)\n', ...
    system_info.total_slosh_mass, 100*system_info.total_slosh_mass/system_info.spacecraft_mass);

% Show tank information
fprintf('\nTank Properties:\n');
for i = 1:system_info.n_tanks
    tank = system_info.tank_info(i);
    fprintf('  %s: mass=%.0fkg, freq=%.2fHz, rate_sens_x=%.2e, rate_sens_z=%.2e\n', ...
        tank.name, tank.m, tank.omega_sl/(2*pi), tank.mu_P_rate_x, tank.mu_P_rate_z);
end

% Compare diagonal vs proper inverse for demonstration
fprintf('\nInertia Matrix Analysis:\n');
fprintf('Original inertia diagonal: [%.0f, %.0f, %.0f]\n', ...
    spacecraft_params.I(1,1), spacecraft_params.I(2,2), spacecraft_params.I(3,3));
fprintf('1/I_diagonal:              [%.2e, %.2e, %.2e]\n', ...
    1/spacecraft_params.I(1,1), 1/spacecraft_params.I(2,2), 1/spacecraft_params.I(3,3));
fprintf('inv(I) diagonal:           [%.2e, %.2e, %.2e]\n', ...
    system_info.inverse_inertia(1,1), system_info.inverse_inertia(2,2), system_info.inverse_inertia(3,3));

if ~system_info.is_diagonal_inertia
    fprintf('Difference is significant due to off-diagonal terms!\n');
else
    fprintf('Values are identical (diagonal inertia matrix)\n');
end

% Check eigenvalues for stability
eigenvals = eig(A);
unstable_modes = find(real(eigenvals) > 1e-6);

if isempty(unstable_modes)
    fprintf('\nSystem is stable - all eigenvalues have negative real parts\n');
else
    fprintf('\nWarning: System has %d unstable mode(s)\n', length(unstable_modes));
end

% Create system object and save to workspace
sys = ss(A, B, C, D);
sys.StateName = state_names;
sys.InputName = input_names;
sys.OutputName = state_names;

assignin('base', 'multi_tank_system', sys);
assignin('base', 'system_info', system_info);
assignin('base', 'A_matrix', A);
assignin('base', 'B_matrix', B);

fprintf('\nSystem saved to workspace as "multi_tank_system"\n');
```

end

%% Analysis function for multi-tank systems
function analyze_multi_tank_system(A, system_info)
fprintf(’\n=== Multi-Tank Slosh System Analysis ===\n’);

```
% Overall system stability
eigenvals = eig(A);
unstable_modes = find(real(eigenvals) > 1e-6);

if isempty(unstable_modes)
    fprintf('✓ Overall system is stable\n');
else
    fprintf('⚠ System has %d unstable mode(s):\n', length(unstable_modes));
    for i = unstable_modes'
        fprintf('  Eigenvalue %d: %.4f + %.4fj\n', i, real(eigenvals(i)), imag(eigenvals(i)));
    end
end

% Individual tank analysis
fprintf('\nIndividual Tank Analysis:\n');
for i = 1:system_info.n_tanks
    tank = system_info.tank_info(i);
    fprintf('  %s:\n', tank.name);
    fprintf('    Natural frequency: %.3f Hz\n', tank.omega_sl/(2*pi));
    fprintf('    Rate sensitivity (x): %.2e rad/s per rad\n', tank.mu_P_rate_x);
    fprintf('    Rate sensitivity (z): %.2e rad/s per rad\n', tank.mu_P_rate_z);
    fprintf('    Slosh mass ratio: %.1f%%\n', 100*tank.m/system_info.spacecraft_mass);
end

% Coupling analysis
fprintf('\nCoupling Strength Analysis:\n');
total_coupling_x = sum([system_info.tank_info.mu_P_rate_x]);
total_coupling_z = sum([system_info.tank_info.mu_P_rate_z]);
fprintf('  Total X-axis coupling: %.2e\n', total_coupling_x);
fprintf('  Total Z-axis coupling: %.2e\n', total_coupling_z);

if max(abs([total_coupling_x, total_coupling_z])) > 0.1
    fprintf('  ⚠ Strong slosh-attitude coupling detected\n');
else
    fprintf('  ✓ Moderate slosh-attitude coupling\n');
end
```

end

% Run example if this file is executed directly
if ~exist(‘spacecraft_params’, ‘var’)
example_multi_tank_usage();
end

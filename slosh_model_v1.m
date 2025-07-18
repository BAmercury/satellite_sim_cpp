function [A, B, C, D, state_names, input_names] = build_slosh_state_space(spacecraft_params, slosh_params)
% BUILD_SLOSH_STATE_SPACE Creates state space model for spacecraft with propellant slosh
%
% Inputs:
%   spacecraft_params: Structure with fields:
%     - Ix, Iy, Iz: Moments of inertia [kg*m^2]
%     - M0: Total spacecraft mass [kg]
%     - g: Effective acceleration (thrust/mass) [m/s^2]
%
%   slosh_params: Structure with fields:
%     - m: Slosh mass [kg]
%     - l: Pendulum length [m]
%     - xH, yH, zH: Hinge location w.r.t. spacecraft CG [m]
%     - zeta: Viscous damping factor [dimensionless]
%
% Outputs:
%   A, B, C, D: State space matrices
%   state_names: Cell array of state variable names
%   input_names: Cell array of input names
%
% State vector: [psi phi theta wx wy wz psi_p theta_p psi_p_dot theta_p_dot]
% Input vector: [Tcx Tcy Tcz fx_thr fy_thr fz_thr]

% Extract spacecraft parameters
Ix = spacecraft_params.Ix;
Iy = spacecraft_params.Iy;
Iz = spacecraft_params.Iz;
M0 = spacecraft_params.M0;
g = spacecraft_params.g;

% Extract slosh parameters
m = slosh_params.m;
l = slosh_params.l;
xH = slosh_params.xH;
yH = slosh_params.yH;
zH = slosh_params.zH;
zeta = slosh_params.zeta;

% Calculate derived parameters
omega_sl = sqrt(g/l * (1 + m/M0));  % Slosh frequency

% Acceleration sensitivities
mu_P_accn_x = m * yH * g / Ix;  % x-axis acceleration sensitivity
mu_P_accn_z = m * yH * g / Iz;  % z-axis acceleration sensitivity

% State vector dimension
N = 10;  % [psi phi theta wx wy wz psi_p theta_p psi_p_dot theta_p_dot]

% Initialize matrices
A = zeros(N, N);
B = zeros(N, 6);
C = eye(N);  % Full state output
D = zeros(N, 6);

%% Build A matrix

% Attitude kinematics: [psi; phi; theta] = [wx; wy; wz]
A(1:3, 4:6) = eye(3);

% Slosh coupling to spacecraft dynamics
% Effect of psi_p on wx (yaw rate)
A(4, 7) = -mu_P_accn_x;
% Effect of theta_p on wz (pitch rate)  
A(6, 8) = -mu_P_accn_z;

% Pendulum kinematics: [psi_p; theta_p] = [psi_p_dot; theta_p_dot]
A(7:8, 9:10) = eye(2);

% Pendulum damping
A(9, 9) = -2*zeta*omega_sl;   % psi_p damping
A(10, 10) = -2*zeta*omega_sl; % theta_p damping

% Pendulum restoring forces (Aψp and Aθp from equations 45-46)
% For psi_p equation (yaw-roll plane)
A(9, 7) = -(1/l) * ((1 + m/M0)*g + (yH - l)*mu_P_accn_z);

% For theta_p equation (roll-pitch plane)  
A(10, 8) = -(1/l) * ((1 + m/M0)*g + (yH - l)*mu_P_accn_x);

%% Build B matrix

% Control torques on spacecraft attitude rates
B(4, 1) = 1/Ix;  % Tcx affects wx
B(5, 2) = 1/Iy;  % Tcy affects wy  
B(6, 3) = 1/Iz;  % Tcz affects wz

% Thruster forces affecting pendulum dynamics
% psi_p equation (from equation 47)
B(9, 2) = -zH/(l*Iy);           % fy_thr effect
B(9, 3) = (yH - l)/(Iz*l);      % fz_thr effect
B(9, 4) = -1/(M0*l);            % fx_thr effect

% theta_p equation  
B(10, 1) = (yH - l)/(Ix*l);     % fx_thr effect
B(10, 2) = -xH/(l*Iy);          % fy_thr effect
B(10, 6) = 1/(M0*l);            % fz_thr effect

%% Define names for clarity
state_names = {‘psi’, ‘phi’, ‘theta’, ‘wx’, ‘wy’, ‘wz’, …
‘psi_p’, ‘theta_p’, ‘psi_p_dot’, ‘theta_p_dot’};

input_names = {‘Tcx’, ‘Tcy’, ‘Tcz’, ‘fx_thr’, ‘fy_thr’, ‘fz_thr’};

end

%% Example usage and validation
function example_usage()
% Define spacecraft parameters
spacecraft_params.Ix = 1000;    % kg*m^2
spacecraft_params.Iy = 2000;    % kg*m^2  
spacecraft_params.Iz = 1500;    % kg*m^2
spacecraft_params.M0 = 1000;    % kg
spacecraft_params.g = 9.81;     % m/s^2 (or effective acceleration)

```
% Define slosh parameters
slosh_params.m = 50;      % kg (slosh mass)
slosh_params.l = 0.5;     % m (pendulum length)
slosh_params.xH = 0.0;    % m (hinge x-location from CG)
slosh_params.yH = 1.0;    % m (hinge y-location from CG)
slosh_params.zH = 0.0;    % m (hinge z-location from CG)
slosh_params.zeta = 0.05; % damping factor

% Build state space model
[A, B, C, D, state_names, input_names] = build_slosh_state_space(spacecraft_params, slosh_params);

% Display results
fprintf('State Space Model Built Successfully!\n');
fprintf('System size: %d states, %d inputs\n', size(A,1), size(B,2));

% Check eigenvalues for stability
eigenvals = eig(A);
fprintf('\nSystem eigenvalues (real parts):\n');
for i = 1:length(eigenvals)
    fprintf('  λ%d = %.4f + %.4fj\n', i, real(eigenvals(i)), imag(eigenvals(i)));
end

% Create system object
sys = ss(A, B, C, D);
sys.StateName = state_names;
sys.InputName = input_names;
sys.OutputName = state_names;

% Plot step response to control torque input
figure;
step(sys(4:6, 1:3), 10);  % Angular rates response to control torques
title('Angular Rate Response to Control Torques');
legend('wx', 'wy', 'wz');
grid on;

% Plot pendulum response to thruster force
figure;
step(sys(7:8, 4:6), 20);  % Pendulum angles response to thruster forces
title('Pendulum Angle Response to Thruster Forces');
legend('ψ_p to f_x', 'θ_p to f_x', 'ψ_p to f_y', 'θ_p to f_y', 'ψ_p to f_z', 'θ_p to f_z');
grid on;

% Return the system for further analysis
assignin('base', 'slosh_system', sys);
assignin('base', 'A_matrix', A);
assignin('base', 'B_matrix', B);
```

end

%% Additional analysis functions

function analyze_slosh_stability(A, spacecraft_params, slosh_params)
% Analyze the stability characteristics of the slosh system

```
eigenvals = eig(A);

% Check for unstable modes
unstable_modes = find(real(eigenvals) > 1e-6);

if isempty(unstable_modes)
    fprintf('System is stable - all eigenvalues have negative real parts\n');
else
    fprintf('Warning: System has %d unstable mode(s)\n', length(unstable_modes));
    for i = unstable_modes'
        fprintf('  Unstable eigenvalue: %.4f + %.4fj\n', real(eigenvals(i)), imag(eigenvals(i)));
    end
end

% Calculate slosh frequency
omega_sl = sqrt(spacecraft_params.g/slosh_params.l * (1 + slosh_params.m/spacecraft_params.M0));
fprintf('\nSlosh natural frequency: %.3f rad/s (%.3f Hz)\n', omega_sl, omega_sl/(2*pi));

% Rate sensitivity
mu_P_rate = (slosh_params.m * slosh_params.yH * spacecraft_params.g / spacecraft_params.Ix) / omega_sl;
fprintf('Rate sensitivity (μ_P_rate): %.6f\n', mu_P_rate);
```

end

function plot_slosh_response(sys, t_final)
% Plot comprehensive slosh response

```
if nargin < 2
    t_final = 30;
end

% Step response to different inputs
figure('Position', [100 100 1200 800]);

% Response to yaw torque
subplot(2,3,1);
step(sys(7:8, 1), t_final);
title('Pendulum Response to Yaw Torque');
ylabel('Angle (rad)');
legend('ψ_p', 'θ_p');
grid on;

% Response to pitch torque
subplot(2,3,2);
step(sys(7:8, 3), t_final);
title('Pendulum Response to Pitch Torque');
ylabel('Angle (rad)');
legend('ψ_p', 'θ_p');
grid on;

% Response to lateral thruster force
subplot(2,3,3);
step(sys(7:8, 4), t_final);
title('Pendulum Response to Lateral Force');
ylabel('Angle (rad)');
legend('ψ_p', 'θ_p');
grid on;

% Angular rate response to pendulum motion (coupling effect)
subplot(2,3,4);
impulse(sys(4, 7), t_final);
title('Yaw Rate Response to ψ_p Impulse');
ylabel('ω_x (rad/s)');
grid on;

subplot(2,3,5);
impulse(sys(6, 8), t_final);
title('Pitch Rate Response to θ_p Impulse');
ylabel('ω_z (rad/s)');
grid on;

% Frequency response
subplot(2,3,6);
bode(sys(7:8, 4:6));
title('Frequency Response: Pendulum/Thruster');
grid on;
```

end

% Run example if this file is executed directly
if ~exist(‘spacecraft_params’, ‘var’)
example_usage();
end

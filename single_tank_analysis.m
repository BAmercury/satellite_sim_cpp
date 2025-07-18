%% OFF-CENTERLINE TANK SLOSH ANALYSIS SCRIPT
% This script demonstrates slosh modeling for a spacecraft with:
% - Single tank offset from vehicle centerline
% - Engines along vehicle centerline (thrust through geometric center)
% - CG may be offset from vehicle centerline due to tank placement
% - Shows impact of off-centerline configuration on slosh coupling

clear; clc; close all;

fprintf(’=== OFF-CENTERLINE TANK SLOSH ANALYSIS ===\n\n’);

%% SPACECRAFT CONFIGURATION
% Define coordinate system:
% X = Roll axis (right wing positive)
% Y = Pitch axis (nose positive, along thrust centerline)  
% Z = Yaw axis (up positive)

% Vehicle geometric centerline (where engines are mounted)
vehicle_centerline = [0, 0, 0];  % Engines thrust through this point

% Spacecraft center of gravity (may be offset due to tank placement)
spacecraft_CG = [0.2, 0.8, 0.1];  % CG is offset from vehicle centerline
fprintf(‘Spacecraft CG offset from vehicle centerline: [%.1f, %.1f, %.1f] m\n’, spacecraft_CG);

% Spacecraft inertia matrix (at the CG, not vehicle centerline!)
% Includes all off-diagonal terms due to asymmetric mass distribution
I_spacecraft = [1200   50   30;    % Asymmetric due to off-centerline tank
50   2000  40;
30    40  1800];

% Total spacecraft mass
M_total = 1500;  % kg

% Effective acceleration (thrust/mass or gravity)
g_eff = 10;  % m/s^2

fprintf(‘Spacecraft inertia matrix (at CG):\n’);
disp(I_spacecraft);
fprintf(‘Total mass: %.0f kg\n’, M_total);

%% TANK CONFIGURATION (Off-Centerline)
% Tank location relative to spacecraft CG
tank_xH = 1.5;   % Tank is 1.5m to the right of CG
tank_yH = 1.0;   % Tank is 1.0m forward of CG
tank_zH = 0.8;   % Tank is 0.8m above CG

% Tank slosh properties
tank_mass = 200;     % kg (sloshing propellant mass)
tank_length = 0.6;   % m (pendulum length)
tank_damping = 0.03; % 3% damping ratio

fprintf(’\nTANK CONFIGURATION:\n’);
fprintf(‘Tank position relative to CG: [%.1f, %.1f, %.1f] m\n’, tank_xH, tank_yH, tank_zH);
fprintf(‘Tank distance from CG: %.2f m\n’, sqrt(tank_xH^2 + tank_yH^2 + tank_zH^2));
fprintf(‘Tank mass (slosh): %.0f kg (%.1f%% of total)\n’, tank_mass, 100*tank_mass/M_total);

% Tank location relative to vehicle centerline (for reference)
tank_position_vehicle = [tank_xH, tank_yH, tank_zH] + spacecraft_CG;
fprintf(‘Tank position relative to vehicle centerline: [%.1f, %.1f, %.1f] m\n’, tank_position_vehicle);

%% ENGINE/THRUSTER CONFIGURATION
% Main engines along vehicle centerline
main_engine_position = [0, 0, 0];  % At vehicle centerline
main_engine_direction = [0, 1, 0]; % Thrust along +Y (pitch axis)

% RCS thrusters for attitude control (example: 4 thrusters in X-Z plane)
rcs_positions = [1.0,  2.0,  1.0;    % Thruster 1: +X, forward, +Z
-1.0,  2.0,  1.0;    % Thruster 2: -X, forward, +Z  
1.0,  2.0, -1.0;    % Thruster 3: +X, forward, -Z
-1.0,  2.0, -1.0];   % Thruster 4: -X, forward, -Z

rcs_directions = [1, 0, 0;     % Thruster 1: +X direction
-1, 0, 0;     % Thruster 2: -X direction
0, 0, 1;     % Thruster 3: +Z direction  
0, 0, -1];   % Thruster 4: -Z direction

rcs_max_force = 50;  % N per thruster

fprintf(’\nENGINE CONFIGURATION:\n’);
fprintf(‘Main engine position (vehicle centerline): [%.1f, %.1f, %.1f] m\n’, main_engine_position);
fprintf(‘Main engine thrust direction: [%.0f, %.0f, %.0f]\n’, main_engine_direction);
fprintf(‘RCS thruster max force: %.0f N each\n’, rcs_max_force);

%% BUILD SPACECRAFT PARAMETERS STRUCTURE
spacecraft_params.I = I_spacecraft;
spacecraft_params.M0 = M_total;
spacecraft_params.g = g_eff;

%% BUILD TANK PARAMETERS STRUCTURE  
tank_params.m = tank_mass;
tank_params.l = tank_length;
tank_params.xH = tank_xH;
tank_params.yH = tank_yH;
tank_params.zH = tank_zH;
tank_params.zeta = tank_damping;
tank_params.name = ‘OffCenterTank’;

%% BUILD SLOSH MODEL
fprintf(’\nBuilding slosh model…\n’);
[A, B, C, D, state_names, input_names] = build_multi_tank_slosh_state_space(spacecraft_params, tank_params);

% Create system
sys_slosh = ss(A, B, C, D);
sys_slosh.StateName = state_names;
sys_slosh.InputName = input_names;
sys_slosh.OutputName = state_names;

fprintf(‘Slosh model built: %d states, %d inputs\n’, size(A,1), size(B,2));

%% CALCULATE COUPLING EFFECTS
fprintf(’\nANALYZING COUPLING EFFECTS:\n’);

% Calculate slosh natural frequency
omega_slosh = sqrt(g_eff/tank_length * (1 + tank_mass/M_total));
freq_slosh_hz = omega_slosh / (2*pi);
fprintf(‘Slosh natural frequency: %.2f Hz\n’, freq_slosh_hz);

% Calculate acceleration sensitivities (enhanced for 3D)
I_inv = inv(I_spacecraft);

% All coupling terms (not just yH!)
mu_accn_roll = tank_mass * tank_yH * g_eff * I_inv(1,1);   % Tank Y-offset affects roll
mu_accn_pitch = tank_mass * tank_xH * g_eff * I_inv(2,2);  % Tank X-offset affects pitch
mu_accn_yaw_from_x = tank_mass * tank_yH * g_eff * I_inv(3,3);  % Tank Y-offset affects yaw
mu_accn_yaw_from_z = tank_mass * tank_xH * g_eff * I_inv(3,3);  % Tank X-offset affects yaw

fprintf(’\nAcceleration sensitivities:\n’);
fprintf(‘Roll coupling (from Y-offset):   %.2e (rad/s²)/(m/s²)\n’, mu_accn_roll);
fprintf(‘Pitch coupling (from X-offset):  %.2e (rad/s²)/(m/s²)\n’, mu_accn_pitch);  
fprintf(‘Yaw coupling (from Y-offset):    %.2e (rad/s²)/(m/s²)\n’, mu_accn_yaw_from_x);
fprintf(‘Yaw coupling (from X-offset):    %.2e (rad/s²)/(m/s²)\n’, mu_accn_yaw_from_z);

% Rate sensitivities
mu_rate_roll = mu_accn_roll / omega_slosh;
mu_rate_pitch = mu_accn_pitch / omega_slosh;
mu_rate_yaw = mu_accn_yaw_from_x / omega_slosh;

fprintf(’\nRate sensitivities:\n’);
fprintf(‘Roll rate sensitivity:   %.2e (rad/s) per rad of slosh\n’, mu_rate_roll);
fprintf(‘Pitch rate sensitivity:  %.2e (rad/s) per rad of slosh\n’, mu_rate_pitch);
fprintf(‘Yaw rate sensitivity:    %.2e (rad/s) per rad of slosh\n’, mu_rate_yaw);

%% EXTRACT SISO TRANSFER FUNCTIONS
fprintf(’\nExtracting SISO transfer functions…\n’);

% Find state indices
wx_idx = find(strcmp(state_names, ‘wx’));
wy_idx = find(strcmp(state_names, ‘wy’));
wz_idx = find(strcmp(state_names, ‘wz’));

% Find input indices  
fx_idx = find(strcmp(input_names, ‘fx_thr’));
fy_idx = find(strcmp(input_names, ‘fy_thr’));
fz_idx = find(strcmp(input_names, ‘fz_thr’));

% Primary coupling: forces to rates
if ~isempty(fx_idx)
H_fx_to_wx = sys_slosh(wx_idx, fx_idx);  % X-force to roll rate
H_fx_to_wy = sys_slosh(wy_idx, fx_idx);  % X-force to pitch rate (cross-coupling!)
H_fx_to_wz = sys_slosh(wz_idx, fx_idx);  % X-force to yaw rate (cross-coupling!)
end

if ~isempty(fy_idx)
H_fy_to_wx = sys_slosh(wx_idx, fy_idx);  % Y-force to roll rate (cross-coupling!)
H_fy_to_wy = sys_slosh(wy_idx, fy_idx);  % Y-force to pitch rate
H_fy_to_wz = sys_slosh(wz_idx, fy_idx);  % Y-force to yaw rate (cross-coupling!)
end

if ~isempty(fz_idx)
H_fz_to_wx = sys_slosh(wx_idx, fz_idx);  % Z-force to roll rate (cross-coupling!)
H_fz_to_wy = sys_slosh(wy_idx, fz_idx);  % Z-force to pitch rate (cross-coupling!)
H_fz_to_wz = sys_slosh(wz_idx, fz_idx);  % Z-force to yaw rate
end

%% ANALYZE CROSS-COUPLING STRENGTH
fprintf(’\nCROSS-COUPLING ANALYSIS:\n’);

% Test frequency for coupling analysis
test_freq = omega_slosh;  % At slosh frequency

if exist(‘H_fx_to_wx’, ‘var’) && exist(‘H_fx_to_wy’, ‘var’)
[mag_main, ~] = bode(H_fx_to_wx, test_freq);     % Main coupling
[mag_cross, ~] = bode(H_fx_to_wy, test_freq);    % Cross coupling

```
cross_coupling_ratio_db = 20*log10(abs(mag_cross)/mag_main);
fprintf('X-force: Cross-coupling (pitch) vs Main (roll) = %.1f dB\n', cross_coupling_ratio_db);

if cross_coupling_ratio_db > -20
    fprintf('  → STRONG cross-coupling! Consider MIMO control design\n');
elseif cross_coupling_ratio_db > -40
    fprintf('  → Moderate cross-coupling\n');
else
    fprintf('  → Weak cross-coupling\n');
end
```

end

%% VISUALIZE THE CONFIGURATION
fprintf(’\nGenerating visualization…\n’);

figure(‘Name’, ‘Off-Centerline Tank Configuration’, ‘Position’, [100 100 1400 800]);

% 3D configuration plot
subplot(2,3,1);
% Plot spacecraft CG
plot3(spacecraft_CG(1), spacecraft_CG(2), spacecraft_CG(3), ‘ko’, ‘MarkerSize’, 10, ‘LineWidth’, 2);
hold on;

% Plot vehicle centerline
plot3(0, 0, 0, ‘bs’, ‘MarkerSize’, 8, ‘LineWidth’, 2);

% Plot tank location
plot3(tank_xH + spacecraft_CG(1), tank_yH + spacecraft_CG(2), tank_zH + spacecraft_CG(3), …
‘ro’, ‘MarkerSize’, 12, ‘LineWidth’, 2);

% Plot RCS thrusters
for i = 1:size(rcs_positions,1)
plot3(rcs_positions(i,1), rcs_positions(i,2), rcs_positions(i,3), ‘g^’, ‘MarkerSize’, 6);
end

% Add lines showing offsets
plot3([spacecraft_CG(1), tank_xH + spacecraft_CG(1)], …
[spacecraft_CG(2), tank_yH + spacecraft_CG(2)], …
[spacecraft_CG(3), tank_zH + spacecraft_CG(3)], ‘r–’, ‘LineWidth’, 1);

plot3([0, spacecraft_CG(1)], [0, spacecraft_CG(2)], [0, spacecraft_CG(3)], ‘k–’, ‘LineWidth’, 1);

xlabel(‘X (m)’); ylabel(‘Y (m)’); zlabel(‘Z (m)’);
title(‘3D Configuration’);
legend(‘Spacecraft CG’, ‘Vehicle Centerline’, ‘Tank’, ‘RCS Thrusters’, ‘Location’, ‘best’);
grid on; axis equal;

% Bode plots of main coupling
subplot(2,3,2);
if exist(‘H_fx_to_wx’, ‘var’)
bode(H_fx_to_wx);
title(‘Main: X-Force → Roll Rate’);
grid on;
end

subplot(2,3,3);
if exist(‘H_fy_to_wy’, ‘var’)
bode(H_fy_to_wy);
title(‘Main: Y-Force → Pitch Rate’);
grid on;
end

% Cross-coupling plots
subplot(2,3,4);
if exist(‘H_fx_to_wy’, ‘var’)
bode(H_fx_to_wy);
title(‘Cross: X-Force → Pitch Rate’);
grid on;
end

subplot(2,3,5);
if exist(‘H_fy_to_wx’, ‘var’)
bode(H_fy_to_wx);
title(‘Cross: Y-Force → Roll Rate’);
grid on;
end

% Comparison plot
subplot(2,3,6);
if exist(‘H_fx_to_wx’, ‘var’) && exist(‘H_fx_to_wy’, ‘var’)
bode(H_fx_to_wx, H_fx_to_wy);
legend(‘Main (Roll)’, ‘Cross (Pitch)’, ‘Location’, ‘best’);
title(‘Main vs Cross Coupling’);
grid on;
end

%% IMPACT ON CONTROL DESIGN
fprintf(’\nIMPACT ON CONTROL DESIGN:\n’);

% Check if slosh frequency is in typical control bandwidth
if freq_slosh_hz < 5
fprintf(’- Slosh frequency (%.1f Hz) is in control bandwidth\n’, freq_slosh_hz);
fprintf(’- Will directly affect control performance\n’);
fprintf(’- Consider slosh damping (baffles) or notch filters\n’);
else
fprintf(’- Slosh frequency (%.1f Hz) above typical control bandwidth\n’, freq_slosh_hz);
fprintf(’- May cause high-frequency oscillations\n’);
end

% Check coupling strength
max_coupling = max([abs(mu_rate_roll), abs(mu_rate_pitch), abs(mu_rate_yaw)]);
if max_coupling > 0.1
fprintf(’- Strong slosh coupling detected (%.2e)\n’, max_coupling);
fprintf(’- Consider MIMO control approach\n’);
elseif max_coupling > 0.01
fprintf(’- Moderate slosh coupling (%.2e)\n’, max_coupling);
fprintf(’- SISO design may work with cross-axis compensation\n’);
else
fprintf(’- Weak slosh coupling (%.2e)\n’, max_coupling);
fprintf(’- SISO design should be adequate\n’);
end

% Tank placement recommendations
fprintf(’\nTANK PLACEMENT ANALYSIS:\n’);
fprintf(’- Current tank distance from CG: %.2f m\n’, sqrt(tank_xH^2 + tank_yH^2 + tank_zH^2));
fprintf(’- Tank mass fraction: %.1f%%\n’, 100*tank_mass/M_total);

if sqrt(tank_xH^2 + tank_yH^2 + tank_zH^2) > 1.0 && tank_mass/M_total > 0.1
fprintf(’- RECOMMENDATION: Significant slosh effects expected\n’);
fprintf(’- Consider: tank baffles, closer CG placement, or active slosh control\n’);
else
fprintf(’- Tank placement is reasonable for slosh control\n’);
end

%% SAVE RESULTS
fprintf(’\nSaving results to workspace…\n’);

% Save key variables to workspace
assignin(‘base’, ‘spacecraft_params’, spacecraft_params);
assignin(‘base’, ‘tank_params’, tank_params);
assignin(‘base’, ‘slosh_system’, sys_slosh);
assignin(‘base’, ‘slosh_frequency_hz’, freq_slosh_hz);
assignin(‘base’, ‘coupling_analysis’, struct(‘mu_rate_roll’, mu_rate_roll, …
‘mu_rate_pitch’, mu_rate_pitch, ‘mu_rate_yaw’, mu_rate_yaw));

if exist(‘H_fx_to_wx’, ‘var’)
assignin(‘base’, ‘H_fx_to_wx’, H_fx_to_wx);
assignin(‘base’, ‘H_fx_to_wy’, H_fx_to_wy);
end

fprintf(’\n=== ANALYSIS COMPLETE ===\n’);
fprintf(‘Key results saved to workspace:\n’);
fprintf(’- slosh_system: Complete state space model\n’);
fprintf(’- H_fx_to_wx, H_fx_to_wy: Main and cross-coupling transfer functions\n’);
fprintf(’- coupling_analysis: Quantitative coupling strengths\n’);

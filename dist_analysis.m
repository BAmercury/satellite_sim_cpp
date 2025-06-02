%% Spacecraft CG Disturbance Torque Analysis
% This script analyzes worst-case disturbance torques caused by main engines
% due to CG offset from nominal thrust vectors

clear; clc; close all;

%% Define System Parameters
% Engine Configuration (example values - replace with your actual data)
n_engines = 3;

% Engine positions relative to spacecraft geometric center [x, y, z] in meters
% Example: 3 engines arranged in triangular pattern
engine_positions = [
    1.5, 0, -2.0;      % Engine 1
    -0.75, 1.3, -2.0;  % Engine 2
    -0.75, -1.3, -2.0  % Engine 3
];

% Engine thrust vectors (unit vectors pointing in thrust direction)
% Assuming all engines thrust in +Z direction (adjust as needed)
engine_thrust_vectors = [
    0, 0, 1;  % Engine 1
    0, 0, 1;  % Engine 2
    0, 0, 1   % Engine 3
];

% Nominal thrust magnitude for each engine (N)
engine_thrust_nominal = [5000; 5000; 5000];  % Equal thrust example

% CG bounds relative to geometric center (meters)
cg_bounds = struct();
cg_bounds.x_min = -0.1;  % Min X offset
cg_bounds.x_max = 0.1;   % Max X offset
cg_bounds.y_min = -0.1;  % Min Y offset
cg_bounds.y_max = 0.1;   % Max Y offset
cg_bounds.z_min = -0.05; % Min Z offset
cg_bounds.z_max = 0.05;  % Max Z offset

% Grid resolution for CG search
n_points_per_axis = 20;  % Increase for finer resolution

%% Create CG Grid
x_cg = linspace(cg_bounds.x_min, cg_bounds.x_max, n_points_per_axis);
y_cg = linspace(cg_bounds.y_min, cg_bounds.y_max, n_points_per_axis);
z_cg = linspace(cg_bounds.z_min, cg_bounds.z_max, n_points_per_axis);

[X_cg, Y_cg, Z_cg] = meshgrid(x_cg, y_cg, z_cg);

% Flatten for easier iteration
cg_positions = [X_cg(:), Y_cg(:), Z_cg(:)];
n_cg_positions = size(cg_positions, 1);

%% Initialize Results Storage
torque_magnitudes = zeros(n_cg_positions, 1);
torque_vectors = zeros(n_cg_positions, 3);

%% Calculate Disturbance Torques for Each CG Position
fprintf('Analyzing %d CG positions...\n', n_cg_positions);

for i = 1:n_cg_positions
    % Current CG position
    cg_current = cg_positions(i, :);
    
    % Calculate total torque for this CG position
    total_torque = [0; 0; 0];
    
    for j = 1:n_engines
        % Moment arm from CG to engine
        r_cg_to_engine = engine_positions(j, :)' - cg_current';
        
        % Thrust force vector
        F_thrust = engine_thrust_nominal(j) * engine_thrust_vectors(j, :)';
        
        % Torque = r × F
        torque_j = cross(r_cg_to_engine, F_thrust);
        total_torque = total_torque + torque_j;
    end
    
    % Store results
    torque_vectors(i, :) = total_torque';
    torque_magnitudes(i) = norm(total_torque);
end

%% Find Worst Case Scenarios
[max_torque_mag, max_idx] = max(torque_magnitudes);
worst_cg = cg_positions(max_idx, :);
worst_torque_vector = torque_vectors(max_idx, :);

% Also find worst case for each axis
[max_torque_x, idx_x] = max(abs(torque_vectors(:, 1)));
[max_torque_y, idx_y] = max(abs(torque_vectors(:, 2)));
[max_torque_z, idx_z] = max(abs(torque_vectors(:, 3)));

%% Display Results
fprintf('\n========== WORST CASE ANALYSIS RESULTS ==========\n');
fprintf('Maximum Total Torque Magnitude: %.2f N-m\n', max_torque_mag);
fprintf('Worst Case CG Position: [%.4f, %.4f, %.4f] m\n', worst_cg);
fprintf('Worst Case Torque Vector: [%.2f, %.2f, %.2f] N-m\n', worst_torque_vector);
fprintf('\nWorst Case Per Axis:\n');
fprintf('  Max |Tx|: %.2f N-m at CG = [%.4f, %.4f, %.4f] m\n', ...
    max_torque_x, cg_positions(idx_x, :));
fprintf('  Max |Ty|: %.2f N-m at CG = [%.4f, %.4f, %.4f] m\n', ...
    max_torque_y, cg_positions(idx_y, :));
fprintf('  Max |Tz|: %.2f N-m at CG = [%.4f, %.4f, %.4f] m\n', ...
    max_torque_z, cg_positions(idx_z, :));

%% Visualization
% 1. 3D scatter plot of torque magnitude vs CG position
figure('Name', 'Torque Magnitude vs CG Position');
scatter3(cg_positions(:,1), cg_positions(:,2), cg_positions(:,3), ...
    50, torque_magnitudes, 'filled');
colorbar;
xlabel('CG X offset (m)');
ylabel('CG Y offset (m)');
zlabel('CG Z offset (m)');
title('Disturbance Torque Magnitude Distribution');
hold on;
% Mark worst case
plot3(worst_cg(1), worst_cg(2), worst_cg(3), 'r*', 'MarkerSize', 20);
text(worst_cg(1), worst_cg(2), worst_cg(3), '  Worst Case', 'Color', 'red');
grid on;

% 2. Slice plots at different Z levels
figure('Name', 'Torque Magnitude Slices');
z_slice_indices = round(linspace(1, n_points_per_axis, 4));
for i = 1:4
    subplot(2,2,i);
    z_idx = z_slice_indices(i);
    torque_slice = reshape(torque_magnitudes, [n_points_per_axis, n_points_per_axis, n_points_per_axis]);
    imagesc(x_cg, y_cg, torque_slice(:,:,z_idx)');
    colorbar;
    xlabel('CG X offset (m)');
    ylabel('CG Y offset (m)');
    title(sprintf('Torque Magnitude at Z = %.3f m', z_cg(z_idx)));
    axis equal tight;
end

% 3. Torque components visualization
figure('Name', 'Torque Components');
subplot(2,2,1);
scatter3(cg_positions(:,1), cg_positions(:,2), cg_positions(:,3), ...
    30, abs(torque_vectors(:,1)), 'filled');
colorbar;
title('|Tx| Distribution');
xlabel('X'); ylabel('Y'); zlabel('Z');

subplot(2,2,2);
scatter3(cg_positions(:,1), cg_positions(:,2), cg_positions(:,3), ...
    30, abs(torque_vectors(:,2)), 'filled');
colorbar;
title('|Ty| Distribution');
xlabel('X'); ylabel('Y'); zlabel('Z');

subplot(2,2,3);
scatter3(cg_positions(:,1), cg_positions(:,2), cg_positions(:,3), ...
    30, abs(torque_vectors(:,3)), 'filled');
colorbar;
title('|Tz| Distribution');
xlabel('X'); ylabel('Y'); zlabel('Z');

%% Monte Carlo Analysis (Optional)
% Run Monte Carlo simulation for more thorough analysis
fprintf('\n\nRunning Monte Carlo analysis...\n');
n_monte_carlo = 10000;
mc_cg_positions = zeros(n_monte_carlo, 3);
mc_torques = zeros(n_monte_carlo, 3);

for i = 1:n_monte_carlo
    % Random CG within bounds
    mc_cg = [
        cg_bounds.x_min + (cg_bounds.x_max - cg_bounds.x_min) * rand();
        cg_bounds.y_min + (cg_bounds.y_max - cg_bounds.y_min) * rand();
        cg_bounds.z_min + (cg_bounds.z_max - cg_bounds.z_min) * rand()
    ];
    
    % Calculate torque
    total_torque = [0; 0; 0];
    for j = 1:n_engines
        r_cg_to_engine = engine_positions(j, :)' - mc_cg;
        F_thrust = engine_thrust_nominal(j) * engine_thrust_vectors(j, :)';
        total_torque = total_torque + cross(r_cg_to_engine, F_thrust);
    end
    
    mc_cg_positions(i, :) = mc_cg';
    mc_torques(i, :) = total_torque';
end

mc_torque_mags = vecnorm(mc_torques, 2, 2);
[mc_max_torque, mc_max_idx] = max(mc_torque_mags);

fprintf('Monte Carlo Results:\n');
fprintf('  Max torque: %.2f N-m\n', mc_max_torque);
fprintf('  At CG: [%.4f, %.4f, %.4f] m\n', mc_cg_positions(mc_max_idx, :));
fprintf('  Mean torque magnitude: %.2f N-m\n', mean(mc_torque_mags));
fprintf('  Std deviation: %.2f N-m\n', std(mc_torque_mags));

% Histogram of torque magnitudes
figure('Name', 'Monte Carlo Torque Distribution');
histogram(mc_torque_mags, 50);
xlabel('Torque Magnitude (N-m)');
ylabel('Frequency');
title('Distribution of Disturbance Torque Magnitudes');
grid on;

%% Export Results
% Save worst case scenarios to file
results = struct();
results.worst_case_cg = worst_cg;
results.worst_case_torque = worst_torque_vector;
results.worst_case_magnitude = max_torque_mag;
results.all_cg_positions = cg_positions;
results.all_torques = torque_vectors;
results.all_magnitudes = torque_magnitudes;
results.engine_config = struct('positions', engine_positions, ...
    'thrust_vectors', engine_thrust_vectors, ...
    'thrust_nominal', engine_thrust_nominal);

save('cg_torque_analysis_results.mat', 'results');
fprintf('\nResults saved to cg_torque_analysis_results.mat\n');

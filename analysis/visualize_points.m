%% Spacecraft Visualization Script
% This script creates a 3D visualization of a spacecraft with:
% - Center of Gravity (CG)
% - Thruster locations
% - Body and ECI reference frames
% - Quaternion-based orientation

clear; close all; clc;

%% Define Spacecraft Properties
% Center of Gravity (in body frame)
CG = [0; 0; 0]; % at origin for simplicity

% Thruster locations (in body frame)
% Example: 8 thrusters at corners of a cube
thruster_positions = [
    1,  1,  1;   % Thruster 1
   -1,  1,  1;   % Thruster 2
   -1, -1,  1;   % Thruster 3
    1, -1,  1;   % Thruster 4
    1,  1, -1;   % Thruster 5
   -1,  1, -1;   % Thruster 6
   -1, -1, -1;   % Thruster 7
    1, -1, -1;   % Thruster 8
]' * 0.5; % Scale to 0.5m

% Spacecraft body outline (simple box)
body_vertices = [
    1,  1,  1;
   -1,  1,  1;
   -1, -1,  1;
    1, -1,  1;
    1,  1, -1;
   -1,  1, -1;
   -1, -1, -1;
    1, -1, -1;
]' * 0.6; % Slightly larger than thrusters

% Define faces for the box
body_faces = [
    1 2 3 4;   % Top
    5 6 7 8;   % Bottom
    1 2 6 5;   % Front
    4 3 7 8;   % Back
    1 4 8 5;   % Right
    2 3 7 6;   % Left
];

%% Initialize Quaternion
% Example quaternion for ECI to body transformation
% Identity quaternion (no rotation)
q_eciTobody = [1; 0; 0; 0]; % [q0; q1; q2; q3]

% You can modify this to test different orientations
% Example: 45 degree rotation about z-axis
% angle = pi/4;
% q_eciTobody = [cos(angle/2); 0; 0; sin(angle/2)];

%% Create Figure
figure('Name', 'Spacecraft Visualization', 'Position', [100 100 800 600]);
hold on; grid on; axis equal;
xlabel('X (m)'); ylabel('Y (m)'); zlabel('Z (m)');
title('Spacecraft Orientation Visualization');
view(30, 20);

%% Function to convert quaternion to rotation matrix
quat2dcm = @(q) [
    q(1)^2+q(2)^2-q(3)^2-q(4)^2,  2*(q(2)*q(3)-q(1)*q(4)),      2*(q(2)*q(4)+q(1)*q(3));
    2*(q(2)*q(3)+q(1)*q(4)),      q(1)^2-q(2)^2+q(3)^2-q(4)^2,  2*(q(3)*q(4)-q(1)*q(2));
    2*(q(2)*q(4)-q(1)*q(3)),      2*(q(3)*q(4)+q(1)*q(2)),      q(1)^2-q(2)^2-q(3)^2+q(4)^2
];

%% Plot ECI Reference Frame (fixed)
scale = 1.5;
quiver3(0, 0, 0, scale, 0, 0, 'r', 'LineWidth', 2);
quiver3(0, 0, 0, 0, scale, 0, 'g', 'LineWidth', 2);
quiver3(0, 0, 0, 0, 0, scale, 'b', 'LineWidth', 2);
text(scale*1.1, 0, 0, 'X_{ECI}', 'Color', 'r', 'FontWeight', 'bold');
text(0, scale*1.1, 0, 'Y_{ECI}', 'Color', 'g', 'FontWeight', 'bold');
text(0, 0, scale*1.1, 'Z_{ECI}', 'Color', 'b', 'FontWeight', 'bold');

%% Apply Rotation to Spacecraft
% Get rotation matrix from quaternion
R_bodyToEci = quat2dcm(q_eciTobody)'; % Transpose for body to ECI

% Transform all spacecraft points
CG_eci = R_bodyToEci * CG;
thruster_positions_eci = R_bodyToEci * thruster_positions;
body_vertices_eci = R_bodyToEci * body_vertices;

%% Plot Spacecraft Body
% Plot body outline
patch('Vertices', body_vertices_eci', 'Faces', body_faces, ...
      'FaceColor', [0.7 0.7 0.7], 'FaceAlpha', 0.3, ...
      'EdgeColor', 'k', 'LineWidth', 1);

%% Plot Center of Gravity
plot3(CG_eci(1), CG_eci(2), CG_eci(3), 'ko', ...
      'MarkerSize', 10, 'MarkerFaceColor', 'k');
text(CG_eci(1)+0.1, CG_eci(2)+0.1, CG_eci(3)+0.1, 'CG', ...
     'FontWeight', 'bold', 'FontSize', 10);

%% Plot Thrusters
for i = 1:size(thruster_positions_eci, 2)
    plot3(thruster_positions_eci(1,i), thruster_positions_eci(2,i), ...
          thruster_positions_eci(3,i), 'rs', 'MarkerSize', 8, ...
          'MarkerFaceColor', 'r');
    text(thruster_positions_eci(1,i)+0.05, ...
         thruster_positions_eci(2,i)+0.05, ...
         thruster_positions_eci(3,i)+0.05, ...
         ['T' num2str(i)], 'FontSize', 8);
end

%% Plot Body Reference Frame
% Body frame axes in ECI coordinates
body_scale = 1;
body_x = R_bodyToEci * [body_scale; 0; 0];
body_y = R_bodyToEci * [0; body_scale; 0];
body_z = R_bodyToEci * [0; 0; body_scale];

quiver3(CG_eci(1), CG_eci(2), CG_eci(3), body_x(1), body_x(2), body_x(3), ...
        'r', 'LineWidth', 2, 'LineStyle', '--');
quiver3(CG_eci(1), CG_eci(2), CG_eci(3), body_y(1), body_y(2), body_y(3), ...
        'g', 'LineWidth', 2, 'LineStyle', '--');
quiver3(CG_eci(1), CG_eci(2), CG_eci(3), body_z(1), body_z(2), body_z(3), ...
        'b', 'LineWidth', 2, 'LineStyle', '--');

text(CG_eci(1)+body_x(1)*1.1, CG_eci(2)+body_x(2)*1.1, CG_eci(3)+body_x(3)*1.1, ...
     'X_{body}', 'Color', 'r', 'FontWeight', 'bold');
text(CG_eci(1)+body_y(1)*1.1, CG_eci(2)+body_y(2)*1.1, CG_eci(3)+body_y(3)*1.1, ...
     'Y_{body}', 'Color', 'g', 'FontWeight', 'bold');
text(CG_eci(1)+body_z(1)*1.1, CG_eci(2)+body_z(2)*1.1, CG_eci(3)+body_z(3)*1.1, ...
     'Z_{body}', 'Color', 'b', 'FontWeight', 'bold');

%% Display Quaternion Information
text_str = sprintf('q_{ECI→body} = [%.3f, %.3f, %.3f, %.3f]', ...
                   q_eciTobody(1), q_eciTobody(2), q_eciTobody(3), q_eciTobody(4));
text(0.02, 0.98, text_str, 'Units', 'normalized', ...
     'FontWeight', 'bold', 'BackgroundColor', 'w');

%% Set axis limits
axis_limit = 2;
xlim([-axis_limit axis_limit]);
ylim([-axis_limit axis_limit]);
zlim([-axis_limit axis_limit]);

%% Add legend
legend({'', '', '', 'ECI Frame', '', '', 'Body', 'CG', 'Thrusters', ...
        'Body Frame'}, 'Location', 'best');

%% Animation Function (optional)
% Uncomment the following section to animate rotation

% % Create time vector
% t = linspace(0, 2*pi, 100);
% 
% for i = 1:length(t)
%     % Update quaternion (example: rotation about z-axis)
%     angle = t(i);
%     q_eciTobody = [cos(angle/2); 0; 0; sin(angle/2)];
%     
%     % Clear current plot
%     cla;
%     
%     % Replot everything with new orientation
%     % [Insert all plotting code here]
%     
%     % Update display
%     drawnow;
%     pause(0.05);
% end

%% Utility Functions for Quaternion Operations

% Function to multiply two quaternions
function q = quatmultiply(q1, q2)
    q = [
        q1(1)*q2(1) - q1(2)*q2(2) - q1(3)*q2(3) - q1(4)*q2(4);
        q1(1)*q2(2) + q1(2)*q2(1) + q1(3)*q2(4) - q1(4)*q2(3);
        q1(1)*q2(3) - q1(2)*q2(4) + q1(3)*q2(1) + q1(4)*q2(2);
        q1(1)*q2(4) + q1(2)*q2(3) - q1(3)*q2(2) + q1(4)*q2(1)
    ];
end

% Function to create quaternion from axis-angle
function q = axis_angle_to_quat(axis, angle)
    axis = axis / norm(axis); % Normalize axis
    q = [cos(angle/2); axis * sin(angle/2)];
end

% Function to normalize quaternion
function q = quatnormalize(q)
    q = q / norm(q);
end

%% Example: Different Orientations
% You can test different orientations by modifying q_eciTobody:

% 90 degree yaw (rotation about z-axis)
% q_eciTobody = [cos(pi/4); 0; 0; sin(pi/4)];

% 45 degree pitch (rotation about y-axis)
% q_eciTobody = [cos(pi/8); 0; sin(pi/8); 0];

% 30 degree roll (rotation about x-axis)
% q_eciTobody = [cos(pi/12); sin(pi/12); 0; 0];

% Combined rotation (yaw-pitch-roll)
% q_yaw = [cos(pi/8); 0; 0; sin(pi/8)];
% q_pitch = [cos(pi/12); 0; sin(pi/12); 0];
% q_roll = [cos(pi/16); sin(pi/16); 0; 0];
% q_eciTobody = quatmultiply(quatmultiply(q_yaw, q_pitch), q_roll);

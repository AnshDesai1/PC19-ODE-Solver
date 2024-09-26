clear; close all; clc;
% Lorenz system parameters
sigma = 10;
rho = 28;
beta = 8/3;

% Time span
t0 = 0;
tf = 50;
tol = 1e-5;

% Initial conditions
N = 50;
initial_conditions = 20 * rand(N, 3) - 10;

% Define the Lorenz system
lorenz = @(t, x) [
    sigma * (x(2) - x(1));
    x(1) * (rho - x(3)) - x(2);
    x(1) * x(2) - beta * x(3);
];

% Setup figure
figure;
hold on;
grid on;
xlabel('x(t)', 'Color', 'w');
ylabel('y(t)', 'Color', 'w');
zlabel('z(t)', 'Color', 'w');
title('Lorenz System', 'Color', 'w');
view(3); % Set view to 3D

az = 35;
el = 25;

set(gca, 'Color', 'k'); % Set the axes background color to black
set(gcf, 'Color', 'k'); % Set the figure background color to black
set(gca, 'XColor', 'w', 'YColor', 'w', 'ZColor', 'w');

xlim([-30 30]);
ylim([-50 50]);
zlim([-10 60]);

view(az,el)
    

trajectories = cell(size(initial_conditions, 1), 1);
maxSteps = 0;
for i = 1:size(initial_conditions, 1)
    fprintf('Iteration: %d\n',i)
    y0 = initial_conditions(i, :)';
    [u, t, counter] = pc113(lorenz, t0, y0, tf, tol, 9);
    trajectories{i} = u; % Store the trajectory
    maxSteps = max(maxSteps, size(u, 1)); % Track the maximum number of steps
end

% Animate all trajectories simultaneously
for step = 1:maxSteps
    fprintf('Step %d of %d\n',step,maxSteps)
    for i = 1:size(initial_conditions, 1)
        if step <= size(trajectories{i}, 1)
            % Plot the current position for each trajectory
            plot3(trajectories{i}(1:step, 1), trajectories{i}(1:step, 2), trajectories{i}(1:step, 3), ...
                'Color', 'w', 'LineWidth', 1.5);
            plot3(trajectories{i}(step, 1), trajectories{i}(step, 2), trajectories{i}(step, 3), ...
                '-', 'Color', 'w', 'MarkerFaceColor', 'w');
        end
    end
    az = az + 1;  % Adjust the rotation speed by changing this increment
    view(az, el);
    drawnow;
end

hold off;

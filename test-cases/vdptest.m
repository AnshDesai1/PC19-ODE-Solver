clear; close all; clc

% Set up time parameters
t0 = 0;
tf = 20;
tol = 1e-5;
tol1 = 2e-5 * tol;

% Range of mu values
mus = linspace(0.10, 3, 100);

% Number of initial conditions
num_initial_conditions = 100;

% Range for random initial conditions
initial_range = [-3, 3];

% Generate random initial conditions once
initial_conditions = initial_range(1) + (initial_range(2) - initial_range(1)) * rand(2, num_initial_conditions);

% Set up the plot
figure(1);
%grid on;
title('Solution of van der Pol Equation');
ylabel('y(t)');
xlabel('x(t)');
xlim([-3 3]);
ylim([-5 5]);
hold on;

% Loop over different values of mu
for i = 1:length(mus)
    mu = mus(i);
    cla;
    f = @(t, y) [y(2); mu * (1 - y(1)^2) * y(2) - y(1)];
    
    % Loop over the pre-generated initial conditions
    for j = 1:num_initial_conditions
        % Use the j-th initial condition
        y0 = initial_conditions(:, j);
    
        % Solve the Van der Pol equation using your custom solver
        [u, t, counter] = pc113(f, t0, y0, tf, tol, 9);
    
        % Plot the trajectory
        plot(u(:,1), u(:,2), 'k', 'LineWidth', 0.5);
    end
    
    % Update title with the current value of mu
    title(sprintf('Solution of van der Pol Equation (\\mu = %.2f)', mu));
    
    % Draw the plot for the current value of mu
    drawnow;
end

hold off;
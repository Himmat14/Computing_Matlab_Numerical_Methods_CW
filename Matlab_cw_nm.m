% MATLAB Code for Part I of the Coursework4
% Numerical Methods for Floating Object Oscillation

%% Q1: Analytical Solution and Eigenvalue Analysis
clc; clear; close all;

% Parameters
rho = 500;          % kg/m^3
rho_l = 1000;       % kg/m^3
g = 9.81;           % m/s^2
L = 1;              % m
y0 = 0.1;           % initial displacement (m)
dy0 = 0;            % initial velocity (m/s)

% Analytical parameters
omega_2 = (rho_l * g) / (rho * L);
omega = sqrt(omega_2);

% Analytical solution
syms t
C1 = y0;
C2 = dy0 / omega;
y_analytical = C1 * cos(omega * t) + C2 * sin(omega * t);

time_range = 0:0.01:10; % Time range for plotting
y_analytical_vals = double(subs(y_analytical, t, time_range));

% Plot analytical solution
figure;
plot(time_range, y_analytical_vals, 'LineWidth', 1.5);
xlabel('Time (s)'); ylabel('Displacement y (m)');
title('Analytical Solution: Displacement vs Time');
grid on;

%% Q1c: First-order ODE System
% Convert second-order ODE to first-order system
% Y = [y; dy/dt], Y_dot = [dy/dt; d2y/dt2]
A = [0, 1; -omega_2, 0];

% Eigenvalues of A
eigenvalues = eig(A);
disp('Eigenvalues of matrix A:');
disp(eigenvalues);

%% Q2a: Explicit Euler Method
% Time settings
dt = 0.1;
t_range = 0:dt:10;

% Initial conditions
Y = [y0; dy0];
Y_euler = zeros(2, length(t_range));
Y_euler(:, 1) = Y;

% Euler integration
for i = 2:length(t_range)
    Y_dot = A * Y_euler(:, i-1);
    Y_euler(:, i) = Y_euler(:, i-1) + dt * Y_dot;
end

% Plot Euler solution
figure;
plot(t_range, Y_euler(1, :), 'LineWidth', 1.5);
xlabel('Time (s)'); ylabel('Displacement y (m)');
title('Explicit Euler: Displacement vs Time');
grid on;

%% Q2b: Runge-Kutta 4th Order Method
Y_rk4 = zeros(2, length(t_range));
Y_rk4(:, 1) = Y;

for i = 2:length(t_range)
    k1 = A * Y_rk4(:, i-1);
    k2 = A * (Y_rk4(:, i-1) + dt/2 * k1);
    k3 = A * (Y_rk4(:, i-1) + dt/2 * k2);
    k4 = A * (Y_rk4(:, i-1) + dt * k3);
    Y_rk4(:, i) = Y_rk4(:, i-1) + dt/6 * (k1 + 2*k2 + 2*k3 + k4);
end

% Plot Runge-Kutta solution
figure;
plot(t_range, Y_rk4(1, :), 'LineWidth', 1.5);
xlabel('Time (s)'); ylabel('Displacement y (m)');
title('Runge-Kutta 4th Order: Displacement vs Time');
grid on;

%% Q2c: MATLAB ode45
ode_func = @(t, Y) [Y(2); -omega_2 * Y(1)];
[t_ode45, Y_ode45] = ode45(ode_func, [0 10], Y);

% Plot ode45 solution
figure;
plot(t_ode45, Y_ode45(:, 1), 'LineWidth', 1.5);
xlabel('Time (s)'); ylabel('Displacement y (m)');
title('ode45: Displacement vs Time');
grid on;

%% Q2d: Comparison of Methods
figure;
plot(time_range, y_analytical_vals, 'k', 'LineWidth', 1.5); hold on;
plot(t_range, Y_euler(1, :), 'r--', 'LineWidth', 1.5);
plot(t_range, Y_rk4(1, :), 'b-.', 'LineWidth', 1.5);
plot(t_ode45, Y_ode45(:, 1), 'g:', 'LineWidth', 1.5);
legend('Analytical', 'Euler', 'Runge-Kutta', 'ode45');
xlabel('Time (s)'); ylabel('Displacement y (m)');
title('Comparison of Solutions');
grid on;

%% Q2e: Velocity vs Position
figure;
plot(Y_euler(1, :), Y_euler(2, :), 'r--', 'LineWidth', 1.5); hold on;
plot(Y_rk4(1, :), Y_rk4(2, :), 'b-.', 'LineWidth', 1.5);
plot(Y_ode45(:, 1), Y_ode45(:, 2), 'g:', 'LineWidth', 1.5);
legend('Euler', 'Runge-Kutta', 'ode45');
xlabel('Position y (m)'); ylabel('Velocity v (m/s)');
title('Velocity vs Position');
grid on;

%% Q2f: Solution Comparison
% Reliability and computational complexity
euler_error = abs(y_analytical_vals(1:length(t_range)) - Y_euler(1, :));
rk4_error = abs(y_analytical_vals(1:length(t_range)) - Y_rk4(1, :));
[~, idx] = min(abs(time_range - t_ode45));
ode45_error = abs(y_analytical_vals(1:idx) - Y_ode45(:, 1)');

figure;
plot(t_range, euler_error, 'r--', 'LineWidth', 1.5); hold on;
plot(t_range, rk4_error, 'b-.', 'LineWidth', 1.5);
plot(t_ode45, ode45_error, 'g:', 'LineWidth', 1.5);
legend('Euler Error', 'Runge-Kutta Error', 'ode45 Error',Location='best');
xlabel('Time (s)'); ylabel('Error');
title('Error Comparison Over Time');
grid on;

%% Q3: Numerical Stability
% Experimenting with different step sizes for Euler
step_sizes = [0.5, 0.05, 0.005];
figure;
hold on;
for dt_exp = step_sizes
    t_range_exp = 0:dt_exp:10;
    Y_euler_exp = zeros(2, length(t_range_exp));
    Y_euler_exp(:, 1) = Y;

    for i = 2:length(t_range_exp)
        Y_dot = A * Y_euler_exp(:, i-1);
        Y_euler_exp(:, i) = Y_euler_exp(:, i-1) + dt_exp * Y_dot;
    end
    plot(t_range_exp, Y_euler_exp(1, :), 'LineWidth', 1.5, 'DisplayName', ['dt = ' num2str(dt_exp)]);
end
xlabel('Time (s)'); ylabel('Displacement y (m)');
title('Explicit Euler with Different Step Sizes');
grid on;
legend(Location="best");

%% Q4: Damping Effect
c = 0.2; % Damping coefficient
A_damped = [0, 1; -omega_2, -c];

% Runge-Kutta 4th Order with Damping
Y_rk4_damped = zeros(2, length(t_range));
Y_rk4_damped(:, 1) = Y;

for i = 2:length(t_range)
    k1 = A_damped * Y_rk4_damped(:, i-1);
    k2 = A_damped * (Y_rk4_damped(:, i-1) + dt/2 * k1);
    k3 = A_damped * (Y_rk4_damped(:, i-1) + dt/2 * k2);
    k4 = A_damped * (Y_rk4_damped(:, i-1) + dt * k3);
    Y_rk4_damped(:, i) = Y_rk4_damped(:, i-1) + dt/6 * (k1 + 2*k2 + 2*k3 + k4);
end

figure;
plot(t_range, Y_rk4_damped(1, :), 'LineWidth', 1.5);
xlabel('Time (s)'); ylabel('Displacement y (m)');
title('Damped Runge-Kutta 4th Order: Displacement vs Time');
grid on;


% Rahul Bhadani

% Four car systems: Follow-the-leader Model

% \xdot_i = v_i
% \vdot_i = beta\cfrac{\Delta v}{s_i^2}
% s_i = x_{i-1} - x_i - L
% \Delta v_i = v_i - v_{i-1}

clc;

beta = 150;
L = 4.0;
N = 4; % Four cars (including leader)

try
    data = readtable("speed.txt");
    t_leader = data.Time;
    v_leader = data.speed;
catch
    t_leader = (0:0.01:50.0)';
    v_leader = 20*(1-exp(-t_leader/5)).*(1 + 0.2*sin(0.2*t_leader));
end


f = figure;
plot(t_leader, v_leader, 'LineWidth',2, 'Color','#254422');
xlabel('Time [s]', 'Interpreter','latex', 'FontSize',14);
ylabel('Speed [m/s]', 'Interpreter','latex', 'FontSize',14)
grid on;
title('Leader Vehicle Speed Profile', 'Interpreter','latex', 'FontSize',16);

%%
% Compute leader position
x_leader = zeros(size(t_leader));

% Use trapezoidal integration
for i = 2: length(t_leader)
    dt = t_leader(i) - t_leader(i-1);
    avg_speed = ( v_leader(i-1) + v_leader(i))/2;
    x_leader(i) = x_leader(i-1) + avg_speed*dt;
end

initial_gap = 13;
v0 = 0.0;

% Initial positions of the followers
x0_followers = zeros(1, N - 1);
x0_followers(1) = x_leader(1) - initial_gap; % Position of the first follower

for k = 2:N-1
    x0_followers(k) = x0_followers(k-1) - initial_gap;
end

% Create Initial State Vector
y0 = zeros(1, 2*(N-1));
y0(1:2:end) = x0_followers; % Positions of all follower vehicles
y0(2:2:end) = v0; % Initial speeds of all follower vehicles

% Solving ODE for followers
t_span = [min(t_leader), max(t_leader)];

follower_ode = @(t, y) multi_follower_dynamics(t, y, t_leader, v_leader, x_leader, beta, L, N);

% Define ODE with nested function for interpolation
[t_sol, y_sol] = ode45(follower_ode, t_span, y0);

x1_sol = interp1(t_leader, x_leader, t_sol, 'pchip');
v1_sol = interp1(t_leader, v_leader, t_sol, 'pchip');


% followers_x_sol = zeros(length(t_sol), N-1);
% followers_v_sol = zeros(length(t_sol), N-1);

follower_states = cell(1, N-1);


for k = 1:N-1
    follower_states{k}.x = y_sol(:, 2*k-1);
    follower_states{k}.v = y_sol(:, 2*k);

end

gaps = zeros(length(t_sol), N)-1;
DeltaV = zeros(length(t_sol), N-1);

gaps(:,1 ) = x1_sol -follower_states{1}.x - L;
DeltaV(:,1) = v1_sol - follower_states{1}.v;


for k = 2:N-1
    gaps(:,k) = follower_states{k-1}.x - follower_states{k}.x - L;
    DeltaV(:,k) = follower_states{k-1}.v - follower_states{k}.v;
end

colors = lines(N);  % Distinct colors for each car

f =figure('Position', [100, 100, 1400, 900]);

% Gap vs time
subplot(2,3,1);
hold on;
for k = 1:(N-1)
    plot(t_sol, gaps(:,k), 'LineWidth', 2, 'Color', colors(k+1,:), ...
        'DisplayName', sprintf('Car %d gap', k+1));
end
xlabel('Time [s]', 'Interpreter', 'latex', 'FontSize', 14);
ylabel('Gap [m]', 'Interpreter', 'latex', 'FontSize', 14);
title('Gaps Between Consecutive Cars', 'Interpreter', 'latex', 'FontSize', 16);
legend('show');
grid on;

subplot(2,3,2);
hold on;
plot(t_leader, v_leader, 'LineWidth', 3, 'Color', colors(1,:), ...
    'DisplayName', 'Leader');
for k = 1:(N-1)
    plot(t_sol, follower_states{k}.v, 'LineWidth', 1.5, 'Color', colors(k+1,:), ...
        'DisplayName', sprintf('Car %d', k+1));
end
xlabel('Time [s]', 'Interpreter', 'latex', 'FontSize', 14);
ylabel('Velocity [m/s]', 'Interpreter', 'latex', 'FontSize', 14);
title('Velocity Profiles', 'Interpreter', 'latex', 'FontSize', 16);
legend('show');
grid on;


subplot(2,3,3);
hold on;
accel_leader = gradient(v_leader, t_leader);
plot(t_leader, accel_leader, 'LineWidth', 3, 'Color', colors(1,:), ...
    'DisplayName', 'Leader');
for k = 1:(N-1)
    accel_follower = gradient(follower_states{k}.v, t_sol);
    plot(t_sol, accel_follower, 'LineWidth', 1.5, 'Color', colors(k+1,:), ...
        'DisplayName', sprintf('Car %d', k+1));
end
xlabel('Time [s]', 'Interpreter', 'latex', 'FontSize', 14);
ylabel('Acceleration [$m/s^2$]', 'Interpreter', 'latex', 'FontSize', 14);
title('Acceleration Profiles', 'Interpreter', 'latex', 'FontSize', 16);
legend('show');
grid on;

% Position trajectories
subplot(2,3,4);
hold on;
plot(t_sol, x1_sol, 'LineWidth', 3, 'Color', colors(1,:), ...
    'DisplayName', 'Leader');
for k = 1:(N-1)
    plot(t_sol, follower_states{k}.x, 'LineWidth', 1.5, 'Color', colors(k+1,:), ...
        'DisplayName', sprintf('Car %d', k+1));
end
xlabel('Time [s]', 'Interpreter', 'latex', 'FontSize', 14);
ylabel('Position [m]', 'Interpreter', 'latex', 'FontSize', 14);
title('Position Trajectories', 'Interpreter', 'latex', 'FontSize', 16);
legend('show', 'Location', 'best');
grid on;

subplot(2,3,5);
hold on;
for k = 1:(N-1)
    plot(t_sol, DeltaV(:,k), 'LineWidth', 1.5, 'Color', colors(k+1,:), ...
        'DisplayName', sprintf('Car %d-\\Delta v', k+1));
end
xlabel('Time [s]', 'Interpreter', 'latex', 'FontSize', 14);
ylabel('$\Delta v$ [m/s]', 'Interpreter', 'latex', 'FontSize', 14);
title('Relative Velocities', 'Interpreter', 'latex', 'FontSize', 16);
legend('show');
grid on;

% Gap vs relative velocity (phase plot)
subplot(2,3,6);
hold on;
for k = 1:(N-1)
    %scatter(gaps(:,k), Deltavs(:,k), 10, t_sol, 'filled');
    plot(gaps(:,k), DeltaV(:,k), 'LineWidth', 1.5, 'Color', colors(k+1,:));
    
end

xlabel('Gap [m]', 'Interpreter', 'latex', 'FontSize', 14);
ylabel('$\Delta v$ [m/s]', 'Interpreter', 'latex', 'FontSize', 14);
title('Gap vs. $\Delta v$ (Colored by Time)', 'Interpreter', 'latex', 'FontSize', 16);
grid on;

sgtitle(sprintf('Follow-the-leader Car-following Model: Continuous Solution using ODEs, %d cars', N), 'Interpreter', 'latex', 'FontSize', 18);


function dydt = multi_follower_dynamics(t, y, t_leader, v_leader, x_leader, beta, L, N)
    
    n_followers = N-1; % number of followers
    dydt =zeros(size(y)); % Solution
    
    % States of the follower
    x_followers = y(1:2:end);
    v_followers = y(2:2:end);
    
    v1 = interp1(t_leader, v_leader, t, 'pchip');
    x1 = interp1(t_leader, x_leader, t, 'pchip');
    
    for k = 1:n_followers
        if k == 1
            x_front = x1;
            v_front = v1;
        else
            x_front = x_followers(k-1);
            v_front = v_followers(k-1);
        end
        gap = x_front - x_followers(k) - L;
    
        % Derivative
        dxdt = v_followers(k);
        dvdt = beta*(v_front - v_followers(k))/ ( gap^2);
    
        dydt(2*k-1) = dxdt; % Position derivative
        dydt(2*k) = dvdt; % Velocity derivative
    end
end

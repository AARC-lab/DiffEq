% Rahul Bhadani

% Two car systems: Follow-the-leader Model


% \xdot_i = v_i
% \vdot_i = beta\cfrac{\Delta v}{s_i^2}
% s_i = x_{i-1} - x_i - L
% \Delta v_i = v_{i-1} - v_i 


beta = 150;
L = 4.0;

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

x_follower0 = -13; % Initial position of the follower
v_follower0 = 0.0; % Initial speed of follower

% Initial condition
y0 = [x_follower0, v_follower0];
t_span = [ min(t_leader), max(t_leader)];

follower_ode = @(t, y) follower_dynamics(t, y, t_leader, v_leader, x_leader, beta, L);


% Define ODE with nested function for interpolation
[t_sol, y_sol] = ode45(follower_ode, t_span, y0);

x1_sol = interp1(t_leader, x_leader, t_sol, 'pchip');
v1_sol = interp1(t_leader, v_leader, t_sol, 'pchip');

% Extract the solution
x2_sol = y_sol(:, 1);
v2_sol = y_sol(:, 2);

s = x1_sol - x2_sol - L;
Deltav = v1_sol - v2_sol;

f = figure;
f.Position = [100, 300, 1500, 800];
subplot(2,3,1);
plot(t_sol, s, 'LineWidth', 2, 'Color', '#FF5733');
xlabel('Time [s]', 'Interpreter', 'latex', 'FontSize', 14);
ylabel('Distance Gap [m]', 'Interpreter', 'latex', 'FontSize', 14);
title('Distance Gap Between Leader and Follower', 'Interpreter', 'latex', 'FontSize', 16);
grid on;
subplot(2,3,2);
plot(t_leader, v_leader, 'LineWidth',2, 'Color','#254422');
hold on;
plot(t_sol, v2_sol, 'LineWidth', 2, 'Color', '#34eb77');
xlabel('Time [s]', 'Interpreter', 'latex', 'FontSize', 14);
ylabel('Velocity [m/s]', 'Interpreter', 'latex', 'FontSize', 14);
title('Velocity of Leader and Follower', 'Interpreter', 'latex', 'FontSize', 16);
legend('Leader', 'Follower', 'Interpreter', 'latex', 'FontSize', 12);
grid on;

accel_follow = gradient(v2_sol, t_sol);
accel_leader = gradient(v_leader, t_leader);
subplot(2,3,3);
hold on;
plot(t_leader, accel_leader, 'LineWidth',2, 'Color','#254422');
plot(t_sol, accel_follow, 'LineWidth',2, 'Color','#34eb77');
xlabel('Time [s]', 'Interpreter', 'latex', 'FontSize', 14);
ylabel('Acceleration [$m/s^2$]', 'Interpreter', 'latex', 'FontSize', 14);
title('Acceleration of the Follower', 'Interpreter', 'latex', 'FontSize', 16);
legend('Leader', 'Follower', 'Interpreter', 'latex', 'FontSize', 12);
grid on;

subplot(2,3,4);
plot(t_sol, x1_sol, 'LineWidth', 2, 'Color', '#4286f4');
hold on;
plot(t_sol, x2_sol, 'LineWidth', 2, 'Color', '#34eb77');
xlabel('Time [s]', 'Interpreter', 'latex', 'FontSize', 14);
ylabel('Position [m]', 'Interpreter', 'latex', 'FontSize', 14);
title('Position of Leader and Follower', 'Interpreter', 'latex', 'FontSize', 16);
legend('Leader', 'Follower', 'Interpreter', 'latex', 'FontSize', 12);
grid on;

subplot(2,3,5);
plot(t_sol, Deltav, 'LineWidth', 2, 'Color', '#FF5733');
xlabel('Time [s]', 'Interpreter', 'latex', 'FontSize', 14);
ylabel('Relative velocity $\Delta v$', 'Interpreter', 'latex', 'FontSize', 14);
title('Relative velocity Between Leader and Follower', 'Interpreter', 'latex', 'FontSize', 16);
grid on;

subplot(2,3,6);
plot(s, Deltav, 'LineWidth', 2, 'Color', '#FF5733');
xlabel('Gap (Relative Distance, $s$)', 'Interpreter', 'latex', 'FontSize', 14);
ylabel('Relative velocity $\Delta v$', 'Interpreter', 'latex', 'FontSize', 14);
title('Relative velocity vs Gap', 'Interpreter', 'latex', 'FontSize', 16);
grid on;

sgtitle('Follow-the-leader Car-following Model: Continuous Solution using ODEs', 'Interpreter', 'latex', 'FontSize', 18);

function dydt = follower_dynamics(t, y, t_leader, v_leader, x_leader, beta, L)
    % Interpolate leader data at time t
    v1 = interp1(t_leader, v_leader, t, 'pchip');
    x1 = interp1(t_leader, x_leader, t, 'pchip');

    % Follower states
    x2 = y(1);
    v2 = y(2);

    gap = (x1 -x2 - L);

    % min_gap = 0.1;
    % if gap <= 0.0
    %     gap = min_gap;
    %     warning('Gap <= 0 at t = %.2f. Using min_gap.', t);
    % end

    % differential equation 
    dx2dt = v2;
    dv2dt = beta*(v1-v2)/(gap^2);
    dydt = [dx2dt; dv2dt];
end


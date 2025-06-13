% Rahul Bhadani

% Two car systems: Follow-the-leader Model


% \xdot_i = v_i
% \vdot_i = beta\cfrac{\Delta v}{s_i^2}
% s_i = x_{i-1} - x_i - L
% \Delta v_i = v_i - v_{i-1}
function follow_the_leader_model()

beta = 51.1;
L = 4.0;

try
    data = readtable("speed.txt");
    t_leader = data.Time;
    v_leader = data.speed;
catch
    t_leader = (0:0.1:50.0)';
    v_leader = 20*(1-exp(-t_leader/5)).*(1 + 0.2*sin(0.2*t_leader));
end


f = figure;
plot(t_leader, v_leader, 'LineWidth',2, 'Color','#254422');
xlabel('Time [s]', 'Interpreter','latex', 'FontSize',14);
xlabel('Speed [m/s]', 'Interpreter','latex', 'FontSize',14)
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

x_follower0 = -10; % Initial position of the follower
v_follower0 = 2.0; % Initial speed of follower

% Initial condition
y0 = [x_follower0, v_follower0];
t_span = [ min(t_leader), max(t_leader)];

% Define ODE with nested function for interpolation
[t_sol, y_sol] = ode45(@follower_ode, t_span, y0);

% Extract the solution
x2_sol = y_sol(:, 1);
v2_sol = y_sol(:, 2);

plot(t_leader, v_leader, 'LineWidth',2, 'Color','#254422');
hold on;
plot(t_sol, v2_sol, 'LineWidth', 2, 'Color', '#34eb77');
xlabel('Time [s]', 'Interpreter', 'latex', 'FontSize', 14);
ylabel('Velocity [m/s]', 'Interpreter', 'latex', 'FontSize', 14);
title('Velocity of Leader and Follower', 'Interpreter', 'latex', 'FontSize', 16);
legend('Leader', 'Follower', 'Interpreter', 'latex', 'FontSize', 12);
grid on;

function dydt = follower_ode(t, y)
    % Interpolate leader data at time t
    v1 = interp1(t_leader, v_leader, t, 'pchip');
    x1 = interp1(t_leader, x_leader, t, 'pchip');

    % Follower states
    x2 = y(1);
    v2 = y(2);

    gap = (x1 -x2 - L);

    min_gap = 0.1;
    if gap <= 0.0
        gap = min_gap;
        warning('Gap <= 0 at t = %.2f. Using min_gap.', t);
    end

    % differential equation 
    dx2dt = v2;
    dv2dt = beta*(v2-v1)/(gap^2);
    dydt = [dx2dt; dv2dt];
end

end

%
% % Initial condition of the leader
% x1_0 = 0;
% v1_0 = v_leader(1);
% 
% % Initial condition of the follower% \xdot_i = v_i
% % \vdot_i = beta\cfrac{\Delta v}{s_i^2}
% % s_i = x_{i-1} - x_i - L
% % \Delta v_i = v_i - v_{i-1}
% x2_0 = -20; % start 20 m behind the leader
% v2_0 = 0; % Start from the rest
% 
% % Combined state vector
% y0 = [x1_0; x2_0; v2_0];
% 
% % time-span for integration
% t_span = [t_leader(1), t_leader(end)];
% t_eval = t_leader; % Evaluate the same time as leader data
% 
% % Solve the ODE system
% options =  odeset('RelTol',1e-6, 'AbsTol',1e-8);
% 
% % Now solve
% [t_sol, y_sol] = ode45( @(t, y) car_following_ODE(t, y, t_leader, v_leader, beta, L), ...
%     t_span, y0, options);
% 
% % Extract positions and velocity from the solution
% x1_sol = y_sol(:, 1);
% x2_sol = y_sol(:, 2);
% v2_sol = y_sol(:, 3);
% 
% v1_sol = v_leader;
% 
% % Calculate the disance gap between cars
% s = x1_sol - x2_sol - L;
% accel_follow = gradient(v2_sol, t_sol);
% 
% f = figure;
% f.Position = [100, 100, 1200, 800];
% 
% 
% subplot(2,3,1);
% plot(t_sol, s, 'LineWidth', 2, 'Color', '#FF5733');
% xlabel('Time [s]', 'Interpreter', 'latex', 'FontSize', 14);
% ylabel('Distance Gap [m]', 'Interpreter', 'latex', 'FontSize', 14);
% title('Distance Gap Between Leader and Follower', 'Interpreter', 'latex', 'FontSize', 16);
% grid on;
% subplot(2,3,2);
% plot(t_eval, v1_sol, 'LineWidth', 2, 'Color', '#4286f4');
% hold on;
% plot(t_sol, v2_sol, 'LineWidth', 2, 'Color', '#34eb77');
% xlabel('Time [s]', 'Interpreter', 'latex', 'FontSize', 14);
% ylabel('Velocity [m/s]', 'Interpreter', 'latex', 'FontSize', 14);
% title('Velocity of Leader and Follower', 'Interpreter', 'latex', 'FontSize', 16);
% legend('Leader', 'Follower', 'Interpreter', 'latex', 'FontSize', 12);
% grid on;
% subplot(2,3,3);
% plot(t_sol, accel_follow, 'LineWidth',2, 'Color','#445378');
% xlabel('Time [s]', 'Interpreter', 'latex', 'FontSize', 14);
% ylabel('Acceleration [m/s^2]', 'Interpreter', 'latex', 'FontSize', 14);
% title('Acceleration of the Follower', 'Interpreter', 'latex', 'FontSize', 16);
% grid on;
% subplot(2,3,4);
% plot(t_sol, x1_sol, 'LineWidth', 2, 'Color', '#4286f4');
% hold on;
% plot(t_sol, x2_sol, 'LineWidth', 2, 'Color', '#34eb77');
% xlabel('Time [s]', 'Interpreter', 'latex', 'FontSize', 14);
% ylabel('Position [m/s]', 'Interpreter', 'latex', 'FontSize', 14);
% title('Position of Leader and Follower', 'Interpreter', 'latex', 'FontSize', 16);
% legend('Leader', 'Follower', 'Interpreter', 'latex', 'FontSize', 12);
% grid on;
% 
% sgtitle('Follow-the-leader Car Following Model', 'FontSize', 18);
% 
% 
% % ODE Function
% 
% function dydt = car_following_ODE(t, y, t_leader, v_leader, beta, L)
% % Unpack the state vector
% x1 = y(1); % Position of leader
% x2 = y(2); % Position of follower
% v2 = y(3); % Velocity of follower
% 
% % Interpolate leader's velocity at time t
% v_leader_t = interp1(t_leader, v_leader, t, 'linear', 'extrap');
% 
% % Compute the gap
% % Compute the gap between the leader and follower
% gap = x1 - x2 - L;
% 
% if gap <= 0.1
%     gap = 0.1;
% end
% 
% dx1dt = v_leader_t; % Leader's position derivative
% dx2dt = v2;
% dv2dt = beta*(v2 - v_leader_t)/(gap^2);
% 
% dydt = [dx1dt; dx2dt; dv2dt];
% 
% end
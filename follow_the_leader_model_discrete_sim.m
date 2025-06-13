% Rahul Bhadani
% Two car systems: Follow-the-leader Model Discrete simulation
% \xdot_i = v_i
% \vdot_i = beta\cfrac{\Delta v}{s_i^2}
% s_i = x_{i-1} - x_i - L
% \Delta v_i = v_i - v_{i-1}

beta = 150;  % Increased from 1.0 for better response
L = 4.0;

try
    data = readtable("speed.txt");
    t_leader = data.Time;
    v_leader = data.speed;
catch
    t_leader = (0:0.01:50.0)';  % Changed from 0.0001 to 0.01 for reasonable time step
    v_leader = 20*(1-exp(-t_leader/5)).*(1 + 0.2*sin(0.2*t_leader));
end

f = figure;
plot(t_leader, v_leader, 'LineWidth',2, 'Color','#254422');
xlabel('Time [s]', 'Interpreter','latex', 'FontSize',14);
ylabel('Speed [m/s]', 'Interpreter','latex', 'FontSize',14);  % Fixed: was xlabel
grid on;
title('Leader Vehicle Speed Profile', 'Interpreter','latex', 'FontSize',16);

% Initial conditions
x_leader = zeros(size(v_leader));
x_follower = zeros(size(v_leader));
v_follower = zeros(size(v_leader));
a_follower = zeros(size(v_leader));

% Set initial conditions
x_leader(1) = 0;
x_follower(1) = -13;  % start 20 m behind the leader
v_follower(1) = 0;    % Start from rest

% Main simulation loop
for i = 2:length(t_leader)
    dt = t_leader(i) - t_leader(i-1);
    
    % Update leader position using previous velocity
    x_leader(i) = x_leader(i-1) + v_leader(i-1) * dt;

    % Use trapezoidal integration
    % avg_speed = ( v_leader(i-1) + v_leader(i))/2;
    % x_leader(i) = x_leader(i-1) + avg_speed*dt;
    
    % Calculate gap (using previous positions for stability)
    s = x_leader(i-1) - x_follower(i-1) - L;
    
    % Ensure minimum gap to avoid division by zero or negative gaps
    % min_gap = 0.5;
    % if s <= min_gap
    %     s = min_gap;
    % end
    
    % Calculate follower acceleration (note: Delta v = v_leader - v_follower)
    Delta_v = v_leader(i-1) - v_follower(i-1);
    a_follower(i) = beta * Delta_v / (s^2);
    
    % Update follower velocity
    v_follower(i) = v_follower(i-1) + a_follower(i) * dt;
    
    % % Ensure non-negative velocity
    % if v_follower(i) < 0
    %     v_follower(i) = 0;
    % end
    
    % Update follower position
    x_follower(i) = x_follower(i-1) + v_follower(i-1) * dt;
end

s = x_leader - x_follower;
Deltav = v_leader - v_follower;

f = figure;
f.Position = [100, 300, 1500, 800];
subplot(2,3,1);
plot(t_leader, s, 'LineWidth', 2, 'Color', '#FF5733');
xlabel('Time [s]', 'Interpreter', 'latex', 'FontSize', 14);
ylabel('Distance Gap [m]', 'Interpreter', 'latex', 'FontSize', 14);
title('Distance Gap Between Leader and Follower', 'Interpreter', 'latex', 'FontSize', 16);
grid on;
subplot(2,3,2);
plot(t_leader, v_leader, 'LineWidth',2, 'Color','#254422');
hold on;
plot(t_leader, v_follower, 'LineWidth', 2, 'Color', '#34eb77');
xlabel('Time [s]', 'Interpreter', 'latex', 'FontSize', 14);
ylabel('Velocity [m/s]', 'Interpreter', 'latex', 'FontSize', 14);
title('Velocity of Leader and Follower', 'Interpreter', 'latex', 'FontSize', 16);
legend('Leader', 'Follower', 'Interpreter', 'latex', 'FontSize', 12);
grid on;

accel_follow = gradient(v_follower, t_leader);
subplot(2,3,3);
plot(t_leader, accel_follow, 'LineWidth',2, 'Color','#445378');
xlabel('Time [s]', 'Interpreter', 'latex', 'FontSize', 14);
ylabel('Acceleration [$m/s^2$]', 'Interpreter', 'latex', 'FontSize', 14);
title('Acceleration of the Follower', 'Interpreter', 'latex', 'FontSize', 16);
grid on;

subplot(2,3,4);
plot(t_leader, v_leader, 'LineWidth', 2, 'Color', '#4286f4');
hold on;
plot(t_leader, v_follower, 'LineWidth', 2, 'Color', '#34eb77');
xlabel('Time [s]', 'Interpreter', 'latex', 'FontSize', 14);
ylabel('Position [m]', 'Interpreter', 'latex', 'FontSize', 14);
title('Position of Leader and Follower', 'Interpreter', 'latex', 'FontSize', 16);
legend('Leader', 'Follower', 'Interpreter', 'latex', 'FontSize', 12);
grid on;



subplot(2,3,5);
plot(t_leader, Deltav, 'LineWidth', 2, 'Color', '#FF5733');
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

sgtitle('Follow-the-leader Car-following Model: Discrete Solution', 'Interpreter', 'latex', 'FontSize', 18);
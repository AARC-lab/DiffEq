% Rahul Bhadani

% Two car systems: Follow-the-leader Model Discrete simulation

% \xdot_i = v_i
% \vdot_i = beta\cfrac{\Delta v}{s_i^2}
% s_i = x_{i-1} - x_i - L
% \Delta v_i = v_i - v_{i-1}

beta = 1.0;
L = 4.0;

try
    data = readtable("speed.txt");
    t_leader = data.Time;
    v_leader = data.speed;
catch
    t_leader = (0:0.0001:50.0)';
    v_leader = 20*(1-exp(-t_leader/5)).*(1 + 0.2*sin(0.2*t_leader));
end


f = figure;
plot(t_leader, v_leader, 'LineWidth',2, 'Color','#254422');
xlabel('Time [s]', 'Interpreter','latex', 'FontSize',14);
xlabel('Speed [m/s]', 'Interpreter','latex', 'FontSize',14)
grid on;
title('Leader Vehicle Speed Profile', 'Interpreter','latex', 'FontSize',16);

% Initial condition of the leader
x1_0 = 0;
v1_0 = v_leader(1);

% Initial condition of the follower
x2_0 = -20; % start 20 m behind the leader
v2_0 = 0; % Start from the rest


x_leader = zeros(size(v_leader));
x_follower = zeros(size(v_leader));
v_follower = zeros(size(v_leader));
a_follower = zeros(size(v_leader));


x_leader(1) = 0;

for  i = 2:length(t_leader)
    x_leader(i) = v_leader(i-1)*(t_leader(i) - t_leader(i-1));
    
    s = x_leader(i) - x_follower(i) - L;

    a_follower(i) = beta*(v_follower(i) - v_leader(i))/(  s*s );

    v_follower(i) = v_follower(i-1) + a_follower(i) * (t_leader(i) - t_leader(i-1));
    x_follower(i) = x_follower(i-1) + v_follower(i-1) * (t_leader(i) - t_leader(i-1));

end

% Plot the follower vehicle speed profile
plot(t_leader, v_leader, 'LineWidth',2, 'Color','#254422');
hold on;
plot(t_leader, v_follower, 'LineWidth', 2, 'Color', '#ff5733');

legend('Leader Speed', 'Follower Speed', 'Location', 'Best');
ylabel('Speed [m/s]', 'Interpreter', 'latex', 'FontSize', 14);
% Method 1: Using dsolve for symbolic solution
syms x(t) t
eqn = diff(x,t) == (t-x)/(t+x);
% Solve the differential equation symbolically
sol_symbolic = dsolve(eqn);
disp('Symbolic solution:');
disp(sol_symbolic);

% Method 2: Numerical solution using ode45
% Define the differential equation as a function
f = @(t,x) (t-x)./(t+x);

% Set initial conditions and time span
t0 = 0.1;  % Start slightly away from t=0 to avoid division by zero
tf = 5;    % End time
x0 = 1;    % Initial condition x(t0) = x0

% Solve numerically
[t_num, x_num] = ode45(f, [t0 tf], x0);

% Plot the solution
figure;
plot(t_num, x_num, 'b-', 'LineWidth', 2);
xlabel('t');
ylabel('x(t)');
title('Solution of dx/dt = (t-x)/(t+x)');
grid on;

% Method 3: Multiple solutions with different initial conditions
figure;
hold on;
initial_conditions = [0.5, 1, 2, 3, 6];
colors = ['r', 'g', 'b', 'm', "k"];

for i = 1:length(initial_conditions)
    x0 = initial_conditions(i);
    [t_num, x_num] = ode45(f, [t0 tf], x0);
    plot(t_num, x_num, colors(i), 'LineWidth', 2, ...
         'DisplayName', ['x_0 = ' num2str(x0)]);
end

xlabel('t');
ylabel('x(t)');
title('Solutions with different initial conditions');
legend('show');
grid on;
hold off;

% Method 4: Phase portrait (if you want to visualize the direction field)
figure;
[T, X] = meshgrid(0.1:0.3:5, -2:0.4:4);
DT = ones(size(T));
DX = (T-X)./(T+X);
% Normalize the arrows
N = sqrt(DT.^2 + DX.^2);
DT = DT./N;
DX = DX./N;

quiver(T, X, DT, DX, 0.5, 'k');
hold on;

% Plot solution curves
for i = 1:length(initial_conditions)
    x0 = initial_conditions(i);
    [t_num, x_num] = ode45(f, [0.1 5], x0);
    plot(t_num, x_num, colors(i), 'LineWidth', 2);
end

xlabel('t');
ylabel('x');
title('Direction field and solution curves');
axis([0.1 5 -2 4]);
grid on;
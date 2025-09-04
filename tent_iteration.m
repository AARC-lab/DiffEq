% Tent Map Graphical Iteration
% Parameters
mu = 0.5;
x_0 = 1/4;
N = 200;

% Define the tent function T(x)
T = @(x, mu) mu.*x.*(x < 0.5 & x >= 0) + mu.*(1 - x).*(x >= 0.5 & x <= 1);

% Create x values for plotting the function
x = 0:0.001:1;
y = T(x, mu);

% Plot the tent function
figure;
hold on;
plot(x, y, 'b-', 'LineWidth', 2, 'DisplayName', 'T(x)');

% Plot the diagonal line y = x
plot([0 1], [0 1], 'k--', 'LineWidth', 1, 'DisplayName', 'y = x');

% Initialize for iteration
x_current = x_0;

% Perform graphical iteration
for n = 1:N
    x_next = T(x_current, mu);
    
    % Draw vertical line from (x_current, 0) to (x_current, T(x_current))
    if n == 1
        line([x_current, x_current], [0, x_next], 'Color', 'red', 'LineWidth', 1);
    else
        line([x_current, x_current], [x_current, x_next], 'Color', 'red', 'LineWidth', 1);
    end
    
    % Draw horizontal line from (x_current, T(x_current)) to (x_next, T(x_current))
    line([x_current, x_next], [x_next, x_next], 'Color', 'red', 'LineWidth', 1);
    
    % Update for next iteration
    x_current = x_next;
    
    % Break if we've reached a fixed point or cycle
    % if abs(x_next - x_current) < 1e-10
    %     break;
    % end
end

% Formatting
xlabel('x');
ylabel('T(x)');
title(sprintf('Graphical Iteration of Tent Map (\\mu = %.1f, x_0 = %.2f, N = %d)', mu, x_0, N));
grid on;
axis([0 1 0 1]);
hold off;

%%
% Tent Map Bifurcation Diagram
clear; clc;

% Define the tent function
T = @(x, mu) mu.*x.*(x < 0.5) + mu.*(1 - x).*(x >= 0.5);

% Parameters for bifurcation diagram
mu_min = 0;
mu_max = 3;
mu_resolution = 2000;  % Number of mu values to test
mu_values = linspace(mu_min, mu_max, mu_resolution);

% Parameters for iteration
x0 = 0.3;  % Initial condition (can be any value in [0,1])
n_transient = 2000;  % Number of iterations to skip (let system settle)
n_plot = 40;        % Number of final iterations to plot

% Storage for bifurcation data
mu_plot = [];
x_plot = [];

fprintf('Computing bifurcation diagram...\n');

% Loop through each mu value
for i = 1:length(mu_values)
    mu = mu_values(i);
    
    % Initialize
    x = x0;
    
    % Skip transient behavior
    for n = 1:n_transient
        x = T(x, mu);
        % Keep x in bounds [0,1]
        if x < 0, x = 0; end
        if x > 1, x = 1; end
    end
    
    % Collect attractor points
    for n = 1:n_plot
        x = T(x, mu);
        % Keep x in bounds [0,1]
        if x < 0, x = 0; end
        if x > 1, x = 1; end
        
        % Store data for plotting
        mu_plot = [mu_plot, mu];
        x_plot = [x_plot, x];
    end
    
    % Progress indicator
    if mod(i, 200) == 0
        fprintf('Progress: %.1f%%\n', 100*i/length(mu_values));
    end
end

% Create bifurcation diagram
figure('Position', [100, 100, 800, 600]);
plot(mu_plot, x_plot, '.', 'MarkerSize', 0.5, 'Color', 'blue');
xlabel('\mu (Parameter)', 'FontSize', 12);
ylabel('x (Attractor Values)', 'FontSize', 12);
title('Bifurcation Diagram for Tent Map: x_{n+1} = T(x_n)', 'FontSize', 14);
grid on;
xlim([mu_min, mu_max]);
ylim([0, 1]);

% Add vertical lines at key bifurcation points
hold on;
plot([1, 1], [0, 1], 'r--', 'LineWidth', 2, 'DisplayName', '\mu = 1 (Critical)');
plot([2, 2], [0, 1], 'r--', 'LineWidth', 2, 'DisplayName', '\mu = 2 (Fully Chaotic)');
legend('show', 'Location', 'northeast');

% Add text annotations
text(0.5, 0.9, 'Fixed Point', 'FontSize', 10, 'Color', 'red', 'FontWeight', 'bold');
text(1.5, 0.9, 'Chaotic Region', 'FontSize', 10, 'Color', 'red', 'FontWeight', 'bold');

hold off;

fprintf('Bifurcation diagram completed!\n');

% Additional analysis: Plot some specific cases
figure('Position', [150, 150, 1200, 400]);

% Case 1: mu = 0.5 (converges to 0)
subplot(1,3,1);
mu1 = 0.5;
x_seq1 = zeros(1, 50);
x_seq1(1) = 0.3;
for i = 2:50
    x_seq1(i) = T(x_seq1(i-1), mu1);
end
plot(1:50, x_seq1, 'bo-', 'LineWidth', 1.5);
title(sprintf('\\mu = %.1f (Fixed Point)', mu1));
xlabel('Iteration n');
ylabel('x_n');
grid on;

% Case 2: mu = 1.5 (chaotic)
subplot(1,3,2);
mu2 = 1.5;
x_seq2 = zeros(1, 50);
x_seq2(1) = 0.3;
for i = 2:50
    x_seq2(i) = T(x_seq2(i-1), mu2);
end
plot(1:50, x_seq2, 'ro-', 'LineWidth', 1.5);
title(sprintf('\\mu = %.1f (Chaotic)', mu2));
xlabel('Iteration n');
ylabel('x_n');
grid on;

% Case 3: mu = 2.0 (fully chaotic)
subplot(1,3,3);
mu3 = 2.0;
x_seq3 = zeros(1, 50);
x_seq3(1) = 0.3;
for i = 2:50
    x_seq3(i) = T(x_seq3(i-1), mu3);
end
plot(1:50, x_seq3, 'go-', 'LineWidth', 1.5);
title(sprintf('\\mu = %.1f (Fully Chaotic)', mu3));
xlabel('Iteration n');
ylabel('x_n');
grid on;

sgtitle('Time Series for Different \mu Values', 'FontSize', 14);
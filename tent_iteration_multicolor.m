% Tent Map Bifurcation Diagram
clear; clc;

% Define the tent function
T = @(x, mu) mu.*x.*(x < 0.5) + mu.*(1 - x).*(x >= 0.5);

% Parameters for bifurcation diagram
mu_min = 0;
mu_max = 2;
mu_resolution = 2000;  % Number of mu values to test
mu_values = linspace(mu_min, mu_max, mu_resolution);

% Parameters for iteration
x0 = 0.3;  % Initial condition (can be any value in [0,1])
n_transient = 1000;  % Number of iterations to skip (let system settle)
n_plot = 100;        % Number of final iterations to plot

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

% Create bifurcation diagram with color-coded regions
figure('Position', [100, 100, 800, 600]);
hold on;

% Define color-coded regions
region1_idx = mu_plot < 1;                    % μ < 1 (Fixed point)
region2_idx = (mu_plot >= 1) & (mu_plot < 2); % 1 ≤ μ < 2 (Chaotic)
region3_idx = mu_plot >= 2;                   % μ = 2 (Fully chaotic)

% Plot each region with different colors
plot(mu_plot(region1_idx), x_plot(region1_idx), '.', 'MarkerSize', 0.8, 'Color', 'blue', 'DisplayName', '\mu < 1 (Fixed Point)');
plot(mu_plot(region2_idx), x_plot(region2_idx), '.', 'MarkerSize', 0.5, 'Color', 'red', 'DisplayName', '1 \leq \mu < 2 (Chaotic)');
plot(mu_plot(region3_idx), x_plot(region3_idx), '.', 'MarkerSize', 0.5, 'Color', 'green', 'DisplayName', '\mu = 2 (Fully Chaotic)');

% Add vertical lines at key bifurcation points
plot([1, 1], [0, 1], 'k--', 'LineWidth', 2, 'DisplayName', '\mu = 1 (Critical)');
plot([2, 2], [0, 1], 'k--', 'LineWidth', 2, 'DisplayName', '\mu = 2 (Boundary)');

% Formatting
xlabel('\mu (Parameter)', 'FontSize', 12);
ylabel('x (Attractor Values)', 'FontSize', 12);
title('Color-Coded Bifurcation Diagram for Tent Map: x_{n+1} = T(x_n)', 'FontSize', 14);
grid on;
xlim([mu_min, mu_max]);
ylim([0, 1]);

% Create legend
legend('show', 'Location', 'northeast', 'FontSize', 10);

% Add colored background regions for clarity
alpha_val = 0.1;  % Transparency
fill([0, 1, 1, 0], [0, 0, 1, 1], 'blue', 'FaceAlpha', alpha_val, 'EdgeColor', 'none', 'HandleVisibility', 'off');
fill([1, 2, 2, 1], [0, 0, 1, 1], 'red', 'FaceAlpha', alpha_val, 'EdgeColor', 'none', 'HandleVisibility', 'off');
fill([2, 2, 2, 2], [0, 0, 1, 1], 'green', 'FaceAlpha', alpha_val, 'EdgeColor', 'none', 'HandleVisibility', 'off');

% Add text annotations with matching colors
text(0.5, 0.9, 'Fixed Point Region', 'FontSize', 11, 'Color', 'blue', 'FontWeight', 'bold', 'HorizontalAlignment', 'center');
text(1.5, 0.9, 'Chaotic Region', 'FontSize', 11, 'Color', 'red', 'FontWeight', 'bold', 'HorizontalAlignment', 'center');
text(2.0, 0.85, 'Fully Chaotic', 'FontSize', 10, 'Color', 'green', 'FontWeight', 'bold', 'HorizontalAlignment', 'center', 'Rotation', 90);

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
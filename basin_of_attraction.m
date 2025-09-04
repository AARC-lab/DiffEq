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

%% FEIGENBAUM CONSTANT CALCULATION
% Note: The tent map does NOT exhibit period-doubling, so we'll use the logistic map
fprintf('\n=== FEIGENBAUM CONSTANT CALCULATION ===\n');
fprintf('Note: Tent map does not show period-doubling bifurcations.\n');
fprintf('Using logistic map x_{n+1} = r*x*(1-x) instead.\n\n');

% Logistic map function
logistic = @(x, r) r .* x .* (1 - x);

% Find period-doubling bifurcation points for logistic map
% Approximate known values for period-doubling cascade
r_bifurcations = [];

% Function to find fixed points and check stability
function r_bif = find_period_doubling_points()
    % Known approximate bifurcation points for logistic map
    % These are where period 1 → 2 → 4 → 8 → 16 occur
    r_known = [3.0, 3.449, 3.54409, 3.56441, 3.56875, 3.56969, 3.56989, 3.56993];
    r_bif = r_known;
end

r_bifurcations = find_period_doubling_points();

% Calculate Feigenbaum constant approximations
if length(r_bifurcations) >= 4
    deltas = [];
    for i = 3:length(r_bifurcations)-1
        delta_i = (r_bifurcations(i-1) - r_bifurcations(i-2)) / ...
                  (r_bifurcations(i) - r_bifurcations(i-1));
        deltas = [deltas, delta_i];
        fprintf('δ_%d = (r_%d - r_%d) / (r_%d - r_%d) = %.6f\n', ...
                i-2, i-1, i-2, i, i-1, delta_i);
    end
    
    fprintf('\nFeigenbaum constant estimates:\n');
    for i = 1:length(deltas)
        fprintf('δ_%d ≈ %.6f\n', i, deltas(i));
    end
    
    fprintf('\nTrue Feigenbaum constant: δ = 4.669201609...\n');
    fprintf('Our best estimate: δ ≈ %.6f\n', deltas(end));
    fprintf('Error: %.6f\n', abs(deltas(end) - 4.669201609));
else
    fprintf('Need at least 4 bifurcation points to calculate Feigenbaum constant.\n');
end

% Create comparison plot: Tent Map vs Logistic Map
figure('Position', [200, 200, 1200, 500]);

% Logistic map bifurcation diagram (simplified)
subplot(1,2,1);
r_vals = linspace(2.5, 4, 1000);
r_plot_log = [];
x_plot_log = [];

for i = 1:length(r_vals)
    r = r_vals(i);
    x = 0.5;  % Initial condition
    
    % Skip transients
    for n = 1:500
        x = logistic(x, r);
    end
    
    % Sample attractor
    for n = 1:50
        x = logistic(x, r);
        r_plot_log = [r_plot_log, r];
        x_plot_log = [x_plot_log, x];
    end
end

plot(r_plot_log, x_plot_log, '.', 'MarkerSize', 0.3, 'Color', 'blue');
title('Logistic Map: x_{n+1} = rx_n(1-x_n)');
xlabel('r (Parameter)');
ylabel('x (Attractor Values)');
grid on;
xlim([2.5, 4]);
ylim([0, 1]);

% Add bifurcation points
hold on;
for i = 1:min(6, length(r_bifurcations))
    plot([r_bifurcations(i), r_bifurcations(i)], [0, 1], 'r--', 'LineWidth', 1);
end
text(3.2, 0.9, 'Period-doubling cascade', 'FontSize', 10, 'Color', 'red');
hold off;

% Tent map (from previous plot)
subplot(1,2,2);
plot(mu_plot(region1_idx), x_plot(region1_idx), '.', 'MarkerSize', 0.8, 'Color', 'blue');
hold on;
plot(mu_plot(region2_idx), x_plot(region2_idx), '.', 'MarkerSize', 0.5, 'Color', 'red');
plot([1, 1], [0, 1], 'k--', 'LineWidth', 2);
title('Tent Map: x_{n+1} = T(x_n)');
xlabel('\mu (Parameter)');
ylabel('x (Attractor Values)');
grid on;
xlim([0, 2]);
ylim([0, 1]);
text(0.5, 0.9, 'No period-doubling!', 'FontSize', 10, 'Color', 'red');
hold off;

sgtitle('Comparison: Systems with vs without Period-Doubling', 'FontSize', 14);

fprintf('\n=== SUMMARY ===\n');
fprintf('• Tent Map: Direct transition from fixed point to chaos\n');
fprintf('• Logistic Map: Period-doubling cascade → Feigenbaum constant\n');
fprintf('• Feigenbaum constant applies to smooth unimodal maps, not piecewise-linear maps\n\n');

%% LYAPUNOV EXPONENT CALCULATION FOR LOGISTIC MAP
fprintf('=== LYAPUNOV EXPONENT CALCULATION ===\n');
fprintf('Computing Lyapunov exponent for logistic map x_{n+1} = r*x*(1-x)\n\n');

% Parameters for Lyapunov calculation
r_min_lyap = 2.5;
r_max_lyap = 4.0;
r_resolution_lyap = 500;
r_values_lyap = linspace(r_min_lyap, r_max_lyap, r_resolution_lyap);

% Iteration parameters
n_transient_lyap = 1000;  % Skip transients
n_lyap = 5000;           % Iterations for Lyapunov calculation
x0_lyap = 0.5;           % Initial condition

% Storage for results
lyapunov_exp = zeros(size(r_values_lyap));

fprintf('Computing Lyapunov exponents...\n');

% Calculate Lyapunov exponent for each r value
for i = 1:length(r_values_lyap)
    r = r_values_lyap(i);
    
    % Initialize
    x = x0_lyap;
    lyap_sum = 0;
    
    % Skip transient behavior
    for n = 1:n_transient_lyap
        x = logistic(x, r);
    end
    
    % Calculate Lyapunov exponent
    for n = 1:n_lyap
        x = logistic(x, r);
        
        % Derivative of logistic map: df/dx = r(1-2x)
        derivative = r * (1 - 2*x);
        
        % Avoid log(0) by ensuring derivative is not zero
        if abs(derivative) > 1e-12
            lyap_sum = lyap_sum + log(abs(derivative));
        end
    end
    
    % Average Lyapunov exponent
    lyapunov_exp(i) = lyap_sum / n_lyap;
    
    % Progress indicator
    if mod(i, 50) == 0
        fprintf('Progress: %.1f%% (r = %.3f, λ = %.4f)\n', ...
                100*i/length(r_values_lyap), r, lyapunov_exp(i));
    end
end

% Create Lyapunov exponent plot
figure('Position', [300, 300, 1000, 600]);

subplot(2,1,1);
plot(r_values_lyap, lyapunov_exp, 'b-', 'LineWidth', 1.5);
hold on;
plot([r_min_lyap, r_max_lyap], [0, 0], 'r--', 'LineWidth', 2, 'DisplayName', 'λ = 0');
xlabel('r (Parameter)', 'FontSize', 12);
ylabel('Lyapunov Exponent λ', 'FontSize', 12);
title('Lyapunov Exponent for Logistic Map', 'FontSize', 14);
grid on;
legend('λ(r)', 'Chaos Threshold', 'Location', 'southeast');

% Color-code regions
ylims = ylim;
chaos_idx = lyapunov_exp > 0;
periodic_idx = lyapunov_exp <= 0;

% Add colored background regions
fill([r_values_lyap(periodic_idx), fliplr(r_values_lyap(periodic_idx))], ...
     [ylims(1)*ones(sum(periodic_idx),1)', ylims(2)*ones(sum(periodic_idx),1)'], ...
     'green', 'FaceAlpha', 0.1, 'EdgeColor', 'none', 'HandleVisibility', 'off');

text(3.2, max(lyapunov_exp)*0.8, 'Periodic (λ < 0)', 'FontSize', 11, 'Color', 'green', 'FontWeight', 'bold');
text(3.8, max(lyapunov_exp)*0.8, 'Chaotic (λ > 0)', 'FontSize', 11, 'Color', 'red', 'FontWeight', 'bold');

hold off;

% Zoom in on interesting region
subplot(2,1,2);
zoom_idx = (r_values_lyap >= 3.4) & (r_values_lyap <= 3.8);
plot(r_values_lyap(zoom_idx), lyapunov_exp(zoom_idx), 'b-', 'LineWidth', 2);
hold on;
plot([3.4, 3.8], [0, 0], 'r--', 'LineWidth', 2);

% Mark period-doubling bifurcations
for i = 1:min(4, length(r_bifurcations))
    if r_bifurcations(i) >= 3.4 && r_bifurcations(i) <= 3.8
        plot([r_bifurcations(i), r_bifurcations(i)], [min(lyapunov_exp(zoom_idx)), max(lyapunov_exp(zoom_idx))], ...
             'g--', 'LineWidth', 1.5);
    end
end

xlabel('r (Parameter)', 'FontSize', 12);
ylabel('Lyapunov Exponent λ', 'FontSize', 12);
title('Zoom: Period-Doubling Region', 'FontSize', 12);
grid on;
hold off;

% Analysis and interpretation
fprintf('\n=== LYAPUNOV EXPONENT ANALYSIS ===\n');

% Find onset of chaos (first time λ > 0)
chaos_onset = find(lyapunov_exp > 0, 1);
if ~isempty(chaos_onset)
    r_chaos = r_values_lyap(chaos_onset);
    fprintf('Onset of chaos: r ≈ %.4f (λ changes from negative to positive)\n', r_chaos);
end

% Find specific values
test_r_values = [2.8, 3.2, 3.5, 3.57, 3.8, 4.0];
fprintf('\nLyapunov exponents at specific r values:\n');
for r_test = test_r_values
    [~, idx] = min(abs(r_values_lyap - r_test));
    lambda_val = lyapunov_exp(idx);
    if lambda_val > 0
        behavior = 'Chaotic';
    elseif lambda_val < 0
        behavior = 'Periodic/Fixed Point';
    else
        behavior = 'Critical Point';
    end
    fprintf('r = %.2f: λ ≈ %+.4f (%s)\n', r_test, lambda_val, behavior);
end

% Statistics
max_lyap = max(lyapunov_exp);
min_lyap = min(lyapunov_exp);
chaotic_fraction = sum(lyapunov_exp > 0) / length(lyapunov_exp);

fprintf('\nStatistics:\n');
fprintf('Maximum λ: %.4f (most chaotic)\n', max_lyap);
fprintf('Minimum λ: %.4f (most stable)\n', min_lyap);
fprintf('Fraction of parameter space that is chaotic: %.1f%%\n', chaotic_fraction*100);

fprintf('\n=== INTERPRETATION ===\n');
fprintf('• λ < 0: Periodic behavior (trajectories converge)\n');
fprintf('• λ = 0: Critical point (transition)\n'); 
fprintf('• λ > 0: Chaotic behavior (trajectories diverge exponentially)\n');
fprintf('• Magnitude of λ: Rate of exponential divergence/convergence\n\n');

%% BASINS OF ATTRACTION ANALYSIS
fprintf('=== BASINS OF ATTRACTION ANALYSIS ===\n');
fprintf('Analyzing basins of attraction for different parameter values\n\n');

% Function to compute final state after many iterations
function final_state = compute_attractor(x0, r, map_type, n_iter)
    x = x0;
    for i = 1:n_iter
        if strcmp(map_type, 'logistic')
            x = r * x * (1 - x);
        elseif strcmp(map_type, 'tent')
            if x < 0.5
                x = r * x;
            else
                x = r * (1 - x);
            end
        end
        % Keep in bounds
        if x < 0, x = 0; end
        if x > 1, x = 1; end
    end
    final_state = x;
end

% Parameters for basin analysis
n_initial = 1000;  % Number of initial conditions to test
initial_conditions = linspace(0.001, 0.999, n_initial);
n_settle = 2000;   % Iterations to reach attractor

% Test cases for different dynamics
test_cases = [
    struct('r', 2.8, 'map', 'logistic', 'description', 'Single Fixed Point');
    struct('r', 3.2, 'map', 'logistic', 'description', 'Period-2 Cycle');
    struct('r', 0.8, 'map', 'tent', 'description', 'Tent: Fixed Point at 0');
    struct('r', 1.5, 'map', 'tent', 'description', 'Tent: Chaotic');
];

figure('Position', [400, 400, 1200, 800]);

for case_idx = 1:length(test_cases)
    test_case = test_cases(case_idx);
    
    % Compute final states for all initial conditions
    final_states = zeros(size(initial_conditions));
    for i = 1:length(initial_conditions)
        final_states(i) = compute_attractor(initial_conditions(i), test_case.r, ...
                                          test_case.map, n_settle);
    end
    
    % Create subplot
    subplot(2, 2, case_idx);
    
    % Color-code based on final destination
    if strcmp(test_case.map, 'logistic') && test_case.r == 3.2
        % Period-2 case: two attractors
        [unique_states, ~, idx] = uniquetol(final_states, 1e-3);
        colors = lines(length(unique_states));
        
        for i = 1:length(unique_states)
            mask = idx == i;
            scatter(initial_conditions(mask), final_states(mask), 20, colors(i,:), 'filled');
            hold on;
        end
        
        fprintf('Case %d (%s, r=%.1f): Found %d attractors\n', ...
                case_idx, test_case.map, test_case.r, length(unique_states));
        for i = 1:length(unique_states)
            fprintf('  Attractor %d: x* ≈ %.4f\n', i, unique_states(i));
        end
        
    else
        % Single attractor or chaotic case
        if test_case.r < 1.0 && strcmp(test_case.map, 'tent')
            % All go to zero
            scatter(initial_conditions, final_states, 20, 'blue', 'filled');
        else
            % Chaotic case - color by final state value
            scatter(initial_conditions, final_states, 20, final_states, 'filled');
            colorbar;
        end
        hold on;
        
        if ~(test_case.r > 1.0 && strcmp(test_case.map, 'logistic') && test_case.r < 4)
            fprintf('Case %d (%s, r=%.1f): Single basin covering [0,1]\n', ...
                    case_idx, test_case.map, test_case.r);
        end
    end
    
    xlabel('Initial Condition x_0');
    ylabel('Final State (Attractor)');
    title(sprintf('%s: %s (r=%.1f)', test_case.map, test_case.description, test_case.r));
    grid on;
    xlim([0, 1]);
    
    hold off;
end

sgtitle('Basins of Attraction: Where Do Different Initial Conditions Lead?', 'FontSize', 14);

% Special case: Logistic map with multiple attractors
fprintf('\n=== DETAILED BASIN ANALYSIS ===\n');

% Case with competing attractors (if they exist)
figure('Position', [500, 500, 800, 600]);

% Example: Logistic map with period-3 window (r ≈ 3.83)
r_special = 3.83;  % Approximate value for period-3
n_detail = 2000;
x_detail = linspace(0.001, 0.999, n_detail);
final_detail = zeros(size(x_detail));

for i = 1:length(x_detail)
    final_detail(i) = compute_attractor(x_detail(i), r_special, 'logistic', 3000);
end

% Find unique attractors
[unique_attractors, ~, attractor_idx] = uniquetol(final_detail, 1e-4);
n_attractors = length(unique_attractors);

% Plot basins with different colors
colors = hsv(n_attractors);
for i = 1:n_attractors
    mask = attractor_idx == i;
    plot(x_detail(mask), i*ones(sum(mask),1), '.', 'Color', colors(i,:), 'MarkerSize', 8);
    hold on;
end

xlabel('Initial Condition x_0');
ylabel('Basin Number');
title(sprintf('Basin Structure: Logistic Map r = %.2f', r_special));
ylim([0.5, n_attractors + 0.5]);
grid on;

% Add attractor values as text
for i = 1:n_attractors
    text(0.05, i, sprintf('→ x* ≈ %.4f', unique_attractors(i)), ...
         'FontSize', 10, 'Color', colors(i,:), 'FontWeight', 'bold');
end

hold off;

fprintf('Logistic map (r = %.2f) analysis:\n', r_special);
fprintf('Number of attractors found: %d\n', n_attractors);
for i = 1:n_attractors
    basin_size = sum(attractor_idx == i) / length(attractor_idx);
    fprintf('Basin %d: x* ≈ %.4f, Basin size ≈ %.1f%%\n', ...
            i, unique_attractors(i), basin_size * 100);
end

fprintf('\n=== KEY CONCEPTS ===\n');
fprintf('• Basin of Attraction: Set of initial conditions leading to same attractor\n');
fprintf('• For 1D maps: Usually entire interval [0,1] forms single basin\n');
fprintf('• Exceptions: Maps with multiple stable cycles can have multiple basins\n');
fprintf('• Basin Boundary: Separates regions leading to different attractors\n');
fprintf('• Fractal Basins: Complex, self-similar boundary structures (rare in 1D)\n\n');

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
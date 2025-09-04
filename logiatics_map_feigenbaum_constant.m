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


%%
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
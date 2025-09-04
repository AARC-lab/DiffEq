% x value
x = 0:0.01:1;
mu = 3.0/2.0; % Choose one value between 0 to 2.
T = mu.*x.*(x < 0.5 & x >= 0) + mu.*(1 - x).*(x >= 0.5 & x <= 1);

%%
T = @(x, mu) mu.*x.*(x < 0.5 & x >= 0) + mu.*(1 - x).*(x >= 0.5 & x <= 1);

% Example usage:
x = 0:0.01:1;
result = T(x, mu);
figure;
plot(x, result)

%%



i = 0:0.01:1;
N = 200;
xnew = zeros(1, N);

xnew(1) = 0.25;
for j = 2:N
    xnew(j) = T(xnew(j-1), mu);
end
figure;
plot(0:1:N-1, xnew);

%%
figure;
plot(x, result)
hold on;
for j =1:length(i)-1
    if j == 1
        line([xnew(j), xnew(j)], [0, xnew(j+1)], 'Color', 'red', 'LineWidth', 2);
    else
        line([xnew(j), xnew(j)], [xnew(j), xnew(j+1)], 'Color', 'red', 'LineWidth', 2);
    end
    line([xnew(j), xnew(j+1)], [xnew(j+1), xnew(j+1)], 'Color', 'blue', 'LineWidth', 2);
end
% Plot the diagonal line y = x
plot([0 1], [0 1], 'k--', 'LineWidth', 1, 'DisplayName', 'y = x');
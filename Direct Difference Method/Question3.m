clear;clc;
format long;
h_values = [0.2; 0.1; 0.05; 0.02; 0.01];
results(length(h_values)) = struct('h', [], 'x', [], 'u1', [], 'u2', [], 'u_true', []);
for k = 1:length(h_values)
    h = h_values(k);
    x = 0:h:1;
    N = length(x) - 1;
    A = zeros(N+2, N+1);
    A(1, 1) = 1/h^2;
    for i = 2:N
       A(i, i-1) = -1/h^2;
       A(i, i) = 2/h^2 + x(i);
       A(i, i+1) = -1/h^2;
    end
    A(N+1, N) = -2/h^2;
    A(N+1, N+1) = x(N+1) + 2/h^2 - 4/h;
    A(N+2, N) = -1/h;
    A(N+2, N+1) = 1/h-2;
    b = zeros(N+2, 1);
    b(1) = 3/h^2;
    for i = 2:N
       b(i) = (2*pi^2 + 2*x(i))*cos(pi*x(i)) + x(i)*(x(i)+1);
    end
    b(N+1) = (2*pi^2 + 2*x(N+1))*cos(pi*x(N+1)) + x(N+1)*(x(N+1)+1) + 2/h;
    b(N+2) = 1;
    u1 = A(1:N+1, :) \ b(1:N+1);
    u2 = A([1:N, N+2], :) \ b([1:N, N+2]);
    results(k).x = x;
    results(k).u1 = u1;
    results(k).u2 = u2;
    results(k).u_true = (2*cos(pi.*x)+x+1)';
end
error_u1 = zeros(length(h_values), 1);
error_u2 = zeros(length(h_values), 1);
for k = 1:length(h_values)
    h = h_values(k);
    x = results(k).x;
    N = length(x) - 1;
    error_u1(k) = sqrt(sum(h*((results(k).u1(2:N) - results(k).u_true(2:N)).^2)));
    error_u2(k) = sqrt(sum(h*((results(k).u2(2:N) - results(k).u_true(2:N)).^2)));
end
Alpha1 = (log(error_u1(2:end)) - log(error_u1(1:end-1))) ./ (log(h_values(2:end)) - log(h_values(1:end-1)));
Alpha2 = (log(error_u2(2:end)) - log(error_u2(1:end-1))) ./ (log(h_values(2:end)) - log(h_values(1:end-1)));
figure;
loglog(h_values(1:end-1), Alpha1, '-o', 'DisplayName', '方式1');
hold on;
loglog(h_values(1:end-1), Alpha2, '-x', 'DisplayName', '方式2');
xlabel('h');
ylabel('误差收敛阶');
title('误差关于网格步长的收敛阶变化图像');
grid on;
legend('show');
hold off;


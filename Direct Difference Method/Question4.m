clear;clc;
format long;
h_values = [0.2; 0.1; 0.05; 0.02; 0.01];
condA = zeros(length(h_values),2);
for k = 1:length(h_values)
    h = h_values(k);
    x = 0:h:1;
    N = length(x) - 1;
    A = zeros(N+2, N+1);
    A(1, 1) = 1;
    A(2, 1) = 1/h^2;
    for i = 2:N
       A(i+1, i-1) = -1/h^2;
       A(i+1, i) = 2/h^2 + x(i);
       A(i+1, i+1) = -1/h^2;
    end
    A(N+2, N) = -2/h^2;
    A(N+2, N+1) = x(N+1) + 2/h^2 - 4/h;
    condA(k,1) = cond(A([1, 3:end], :));
    condA(k,2) = cond(A(2:end, :));
end
v1 = -(log(condA(2:end,1)) - log(condA(1:end-1,1))) ./ (log(h_values(2:end)) - log(h_values(1:end-1)));
v2 = -(log(condA(2:end,2)) - log(condA(1:end-1,2))) ./ (log(h_values(2:end)) - log(h_values(1:end-1)));
figure;
loglog(h_values, condA(:,1), '-o', 'DisplayName', '方式1');
hold on;
loglog(h_values, condA(:,2), '-x', 'DisplayName', '方式2');
xlabel('h');
ylabel('离散矩阵A的条件数');
title('条件数关于网格步长的变化图像');
grid on;
legend('show');
figure;
loglog(h_values(1:end-1), v1, '-o', 'DisplayName', '方式1');
hold on;
loglog(h_values(1:end-1), v2, '-x', 'DisplayName', '方式2');
xlabel('h');
ylabel('离散矩阵A的条件数增长速率');
title('条件数关于网格步长的增长速率变化图像');
grid on;
legend('show');
hold off;


clc;
% n = input("输入面元数量：");
n = 100;

sit = zeros(n, 2);
w = zeros(n, 1);
for i=1:n
    w(i) = 2*i*pi/n;
    sit(i,1) = cos(w(i));
    sit(i,2) = sin(w(i));
end

Cp = FindCp(sit, n);

dealt = 0;
for i = 1:n
    dealt = dealt + abs(Cp(i) - (1 - 4*sin(w(i))^2));
end
dealt = dealt / n / 4 * 100;
disp(['平均相对误差为： ', num2str(dealt), ' %']);

% scatter(w,Cp)
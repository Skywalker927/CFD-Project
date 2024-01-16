function Cp = FindCp (sit, n)
Vinf = 1;  %无穷远来流速度

%线源的两端，第n到第1个点为第1个线源，1到2点为第2个线源，以此类推。
xy = zeros(2, 2, n);
for i=1:n
    xy(:, :, i) = [sit(i, :); sit(mod(i, n)+1, :)]';
end

re = zeros(n, n, 2);
a = sqrt(sum((xy(:, 1, :)-xy(:, 2, :)).^2))./2;
for i=1:n
    %i面元中点
    ic = (xy(:, 1, i) + xy(:, 2, i))./2;
    for j=1:n
        if i ~= j
            %j面元中点
            jc = (xy(:, 1, j) + xy(:, 2, j))./2;
            %j面元矢量化，逆时针为正
            fj = xy(:, 2, j) - xy(:, 1, j);
            gg = ic - jc;
            %构造坐标变换过渡矩阵
            fjl = fj./2./a(i);
            fjn = [0, 1; -1, 0] * fjl;
            tranj = [fjl, fjn];
            %j对i的扰动
            gf = tranj \ gg;
            v = [1/4/pi*(log((gf(1)+a(i))^2 + gf(2)^2) - log((gf(1)-a(i))^2 + gf(2)^2));
                1/2/pi*(atan((gf(1)+a(i))/gf(2)) - atan((gf(1)-a(i))/gf(2)))];
            %转为全局坐标
            re(i, j, :) = tranj * v;
        end
    end
end
%j面元对i面元的切向和法向速度
chre = zeros(n, n, 2);
for i=1:n
    fi = xy(:, 2, i) - xy(:, 1, i);
    fil = fi./2./a(i);
    fin = [0, 1; -1, 0] * fil;
    for j=1:n
        if i ~= j
            chre(i, j, :) = [re(i, j, 1), re(i, j, 2)] * [fin, fil];
        end
    end
end
%来流在面元上的切向和法向速度
vre = zeros(n, 2);
for i=1:n
    fi = xy(:, 2, i) - xy(:, 1, i);
    fil = fi./2./a(i);
    fin = [0, 1; -1, 0] * fil; 
    vre(i, :) = [fin(1), fil(1)] .* Vinf;
end
%方程组的系数矩阵
M = chre(:, :, 1) + 0.5 * eye(n);

%面元的强度
x = linsolve(M, -vre(:, 1));

%面元的表面速度
Voo = chre(:, :, 2) * x + vre(:, 2);

Cp = 1 - (Voo./Vinf).^2;

end
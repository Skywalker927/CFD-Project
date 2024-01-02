function Cp = FindCp (sit, n)
Vinf = 1;  %无穷远来流速度

%线源的两端，第n到第1个点为第1个线源，1到2点为第2个线源，以此类推。
xy = zeros(n, 4);
xy(1, 1) = sit(n, 1);
xy(1, 2) = sit(n, 2);
xy(1, 3) = sit(1, 1);
xy(1, 4) = sit(1, 2);
for i=2:n
    xy(i, 1) = sit(i-1, 1);
    xy(i, 2) = sit(i-1, 2);
    xy(i, 3) = sit(i, 1);
    xy(i, 4) = sit(i, 2);
end

re = zeros(n, n, 2);
for i=1:n
    a = sqrt((xy(i, 1)-xy(i, 3))^2 + (xy(i, 2)-xy(i, 4))^2)/2;

    %i面元中点
    ic = zeros(2, 1);
    ic(1) = (xy(i, 1) + xy(i, 3))/2;
    ic(2) = (xy(i, 2) + xy(i, 4))/2;
    for j=1:n
        if i ~= j
            %j面元中点
            jc = zeros(2, 1);
            jc(1) = (xy(j, 1) + xy(j, 3))/2;
            jc(2) = (xy(j, 2) + xy(j, 4))/2;
            %j面元矢量化，逆时针为正
            fj = zeros(2, 1);
            fj(1) = xy(j, 3) - xy(j, 1);
            fj(2) = xy(j, 4) - xy(j, 2);
            %j的中点指向i的中点的矢量
            gg = ic - jc;
            %构造坐标变换过渡矩阵
            fjl = fj'./2./a;
            fjn = (fj' * [0, -1; 1, 0])./2./a;
            tranj = [fjl', fjn'];
            %j对i的扰动
            gf = tranj \ gg;            
            v = zeros(2, 1);
            v(1) = 1/4/pi*(log((gf(1)+a)^2 + gf(2)^2) - log((gf(1)-a)^2 + gf(2)^2));
            v(2) = 1/2/pi*(atan((gf(1)+a)/gf(2)) - atan((gf(1)-a)/gf(2)));            
            %转为全局坐标    
            vv = tranj * v;

            re(i, j, 1) = vv(1);
            re(i, j, 2) = vv(2);
        end
    end
end
%j面元对i面元的切向和法向速度
chre = zeros(n, n, 2);
for i=1:n
    fi = zeros(2, 1);
    fi(1) = xy(i, 3) - xy(i, 1);
    fi(2) = xy(i, 4) - xy(i, 2);
    fil = fi'./2./a;
    fin = (fi' * [0, -1; 1, 0])./2./a;
    for j=1:n
        if i ~= j
            chre(i, j, 1) = fin(2)*re(i, j, 2) + fin(1)*re(i, j, 1);
            chre(i, j, 2) = fil(2)*re(i, j, 2) + fil(1)*re(i, j, 1);
        end
    end
end
%来流在面元上的切向和法向速度
vre = zeros(n, 2);
for i=1:n
    fi = zeros(2, 1);
    fi(1) = xy(i, 3) - xy(i, 1);
    fi(2) = xy(i, 4) - xy(i, 2);
    fil = fi'./2./a;
    fin = (fi' * [0, -1; 1, 0])./2./a; 
    vre(i, 1) = fin(1)*Vinf;
    vre(i, 2) = fil(1)*Vinf;
end
%方程组的系数矩阵
M = zeros(n, n);
for i=1:n
    M(i, i) = 1/2;
    for j=1:n
        if i ~= j
            M(i, j) = chre(i, j, 1);
        end
    end
end

%高斯消元解线性方程组
resu = zeros(n, 1);
for i=1:n
    resu(i) = -1*vre(i, 1);
end

for i=1:n-1
    for j=i+1:n
        bb = M(j, i)/M(i, i);
        for k=i:n
            M(j, k) = M(j, k) - M(i, k)*bb;
        end
        resu(j) = resu(j) - resu(i)*bb;
    end
end
%面元的强度
x = zeros(n, 1);
for i=n:-1:1
    x(i) = resu(i);
    for j=n:-1:i+1
        x(i) = x(i) - M(i, j)*x(j);
    end
    x(i) = x(i)/M(i, i);
end
%面元的表面速度
Voo = zeros(n, 1);
for i=1:n
    Voo(i) = 0;
    for j=1:n
        if i ~= j
            Voo(i) = Voo(i) + x(j)*chre(i, j, 2);
        end
    end
    Voo(i) = Voo(i) + vre(i, 2);
end

Cp = 1 - (Voo./Vinf).^2;

end
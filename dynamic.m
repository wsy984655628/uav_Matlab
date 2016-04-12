function y = dynamic(t, x)

%全局变量
global vin p m;
vin = 10;
p = 1.292;
m = 3;

%变量初始化
pb = x(1);  qb = x(2);  rb = x(3);
ue = x(4);  ve = x(5);  we = x(6);
fi = x(7);  the = x(8); pus = x(9);

ne = 0;
Fi1 = 0;    Fi2 = 0;    Fi3 = 0;    Fi4 = 0;

%参数初始化
Ixx = 0.26;
Iyy = 0.32;
Izz = 0.35;

g = 9.8;

%坐标转换
Tbe = [1 0 0;0 cos(fi) sin(fi);0 -sin(fi) cos(fi);] * ...
        [cos(the) 0 -sin(the);0 1 0;sin(the) 0 cos(the);] * ...
        [cos(pus) sin(pus) 0;-sin(pus) cos(pus) 0;0 0 1;];

Vb = Tbe * [ue ve we]';
ub = Vb(1); vb = Vb(2); wb = Vb(3);

%螺旋桨
R = 0.1905;     %15 inch
A = pi * R^2;
Tp = 3 * quad(@propellert, 0.08, R);
Mp = 3 * quad(@propellerm, 0.08, R);
vin = -wb / 2 + ((wb / 2)^2 + Tp / (2 * p * A)) ^ 0.5;

%舵片
Omg = [pb, qb, rb]';
[Fru, Mru] = rudder(Omg, Vb, Fi1, Fi2, Fi3, Fi4);

% 导流片
Fi = 10;
[Frdi, Mrdi] = diversion(Vb, Fi);

%机身   粗略计算，搭建框架，待详细建模
S = 0.5;
Fbody = -0.5 * p * S * ((wb + vin)^2 + ub^2 + vb^2)^0.5 * [ub vb wb]' * 0.1;

%重力
Fg = Tbe * m * [0 0 g]';

%合力
F = Fru + Frdi + Fbody + Fg + [0 0 Tp]';
M = Mru + Mrdi + [0 0 Mp]';

y(1) = M(1) / Ixx - qb * rb * (Izz - Iyy) / Ixx;
y(2) = M(2) / Iyy - pb * rb * (Ixx - Izz) / Iyy;
y(3) = M(3) / Izz - pb * qb * (Iyy - Ixx) / Izz;

invub = F(1) / m + vb * rb - wb * qb;
invvb = F(2) / m + wb * pb - ub * rb;
invwb = F(3) / m + ub * qb - vb * pb;

y(4:6) = Tbe \ [invub invvb invwb]';

y(7:9) = [1, tan(the)*sin(fi), cos(fi)*tan(the);
          0, cos(fi),          -sin(fi);
          0, sin(fi)/cos(the), cos(fi)/cos(the)] * [pb qb rb]';
end


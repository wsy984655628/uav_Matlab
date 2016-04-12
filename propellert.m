function [ dT ] = propellert( r )
%初始化
global ns Vb vin p
b = 0.04;    %弦长
h = 0.127;  %0.5 inch
wb = Vb(3);

%安装角计算
Rig = atan( h./ (2 * pi * r)) * 180 / pi;

%微元分析
W = ((wb + vin)^2 + (2 * pi * r * ns).^2 ).^0.5;
bet = atan((wb + vin)./(2 * pi * r * ns));
alf = Rig - bet * 180 / pi;

[cl, cd] = CLARK(alf);

dD = 0.5 * b * p * cd.* W.^2;
dL = 0.5 * b * p * cl.* W.^2;

dT = dL.* cos(bet) - dD.* sin(bet);
end


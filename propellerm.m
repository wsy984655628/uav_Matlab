function [ dM ] = propellerm( r )
%��ʼ��
global ns Vb vin p
b = 0.04;    %�ҳ�
h = 0.0127;  %0.5 inch
wb = Vb(3);

%��װ�Ǽ���
Rig = atan( h / (2 * pi * r));

%΢Ԫ����
W = ((wb + vin)^2 + (2 * pi * r * ns).^2 ).^0.5;
bet = atan((wb + vin)./(2 * pi * r * ns));
alf = Rig - bet * 180 / pi;

[cl, cd] = CLARK(alf);

dD = 0.5 * b * p * cd.* W.^2;
dL = 0.5 * b * p * cl.* W.^2;

dM = -r.* (dL.* sin(bet) + dD.* cos(bet));
end


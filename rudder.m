function [ Fru, Mru ] = rudder( Omg, Vb, Fi1, Fi2, Fi3, Fi4 )

%≥ı ºªØ
global vin p;
L = 0.135;
l = 0.12;
S = 0.132 * 0.045;

T1 = [0 L 0; L 0 l; 0 -l 0];
T2 = [0 L -l; L 0 0; l 0 0];
T3 = [0 L 0; L 0 -l; 0 l 0];
T4 = [0 L l; L 0 0; -l 0 0];

V1 = T1 * Omg + Vb;
V2 = T2 * Omg + Vb;
V3 = T3 * Omg + Vb;
V4 = T4 * Omg + Vb;

u1 = V1(1); v1 = V1(2); w1 = V1(3);
u2 = V2(1); v2 = V2(2); w2 = V2(3);
u3 = V3(1); v3 = V3(2); w3 = V3(3);
u4 = V4(1); v4 = V4(2); w4 = V4(3);

bet1 = atan(v1 / (vin + w1));
bet2 = -atan(u2 / (vin + w2));
bet3 = -atan(v3 / (vin + w3));
bet4 = atan(u4 / (vin + w4));

alf1=Fi1+bet1*180/pi;
alf2=Fi2+bet2*180/pi;
alf3=Fi3+bet3*180/pi;
alf4=Fi4+bet4*180/pi;

V1abs = (v1^2 + (w1 + vin)^2)^0.5;
V3abs = (v3^2 + (w3 + vin)^2)^0.5;
V2abs = (u2^2 + (w2 + vin)^2)^0.5;
V4abs = (u4^2 + (w4 + vin)^2)^0.5;

[cl1, cd1] = CLARKy(alf1);
[cl2, cd2] = CLARKy(alf2);
[cl3, cd3] = CLARKy(alf3);
[cl4, cd4] = CLARKy(alf4);

D1 = 0.5 * p * S * V1abs^2 * cd1;
D2 = 0.5 * p * S * V2abs^2 * cd2;
D3 = 0.5 * p * S * V3abs^2 * cd3;
D4 = 0.5 * p * S * V4abs^2 * cd4;

L1 = 0.5 * p * S * V1abs^2 * cl1;
L2 = 0.5 * p * S * V2abs^2 * cl2;
L3 = 0.5 * p * S * V3abs^2 * cl3;
L4 = 0.5 * p * S * V4abs^2 * cl4;

F1 = L1 * cos(bet1) + D1 * sin(bet1);
F2 = L2 * cos(bet2) + D2 * sin(bet2);
F3 = L3 * cos(bet3) + D3 * sin(bet3);
F4 = L4 * cos(bet4) + D4 * sin(bet4);

Q1 = L1 * sin(bet1) - D1 * cos(bet1);
Q2 = L2 * sin(bet2) - D2 * cos(bet2);
Q3 = L3 * sin(bet3) - D3 * cos(bet3);
Q4 = L4 * sin(bet4) - D4 * cos(bet4);

Fru = [(F2 - F4), (F3 - F1), (Q1 + Q2 + Q3 + Q4)]';
Mru = [(L * (F3 - F1) + l * (Q2 - Q4)), (L * (F4 - F2) + l * (Q3 - Q1)), (l*(F1 + F2 + F3 + F4))]';

end


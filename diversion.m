function [ Frdi, Mrdi ] = diversion( Vb, Fi )

%≥ı ºªØ
global vin p;
l = 0.12;
S = 0.147 * 0.06;

V = Vb(3) + vin;

Vabs = abs(V);

alf = Fi;

[cl, cd] = CLARKy(alf);

D = 0.5 * p * S * Vabs^2 * cd;
L = 0.5 * p * S * Vabs^2 * cl;

F = 8 * L;
Q = 8 * D;

Frdi = [0 0 Q]';
Mrdi = [0 0 l * F]';

end


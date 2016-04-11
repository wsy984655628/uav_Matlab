function [sys,x0,str,ts] = hybird(t,x,u,flag)
%SFUNTMPL General MATLAB S-Function Template
%   With MATLAB S-functions, you can define you own ordinary differential
%   equations (ODEs), discrete system equations, and/or just about
%   any type of algorithm to be used within a Simulink block diagram.
%
%   The general form of an MATLAB S-function syntax is:
%       [SYS,X0,STR,TS,SIMSTATECOMPLIANCE] = SFUNC(T,X,U,FLAG,P1,...,Pn)
%
%   What is returned by SFUNC at a given point in time, T, depends on the
%   value of the FLAG, the current state vector, X, and the current
%   input vector, U.
%
%   FLAG   RESULT             DESCRIPTION
%   -----  ------             --------------------------------------------
%   0      [SIZES,X0,STR,TS]  Initialization, return system sizes in SYS,
%                             initial state in X0, state ordering strings
%                             in STR, and sample times in TS.
%   1      DX                 Return continuous state derivatives in SYS.
%   2      DS                 Update discrete states SYS = X(n+1)
%   3      Y                  Return outputs in SYS.
%   4      TNEXT              Return next time hit for variable step sample
%                             time in SYS.
%   5                         Reserved for future (root finding).
%   9      []                 Termination, perform any cleanup SYS=[].
%
%
%   The state vectors, X and X0 consists of continuous states followed
%   by discrete states.
%
%   Optional parameters, P1,...,Pn can be provided to the S-function and
%   used during any FLAG operation.
%
%   When SFUNC is called with FLAG = 0, the following information
%   should be returned:
%
%      SYS(1) = Number of continuous states.
%      SYS(2) = Number of discrete states.
%      SYS(3) = Number of outputs.
%      SYS(4) = Number of inputs.
%               Any of the first four elements in SYS can be specified
%               as -1 indicating that they are dynamically sized. The
%               actual length for all other flags will be equal to the
%               length of the input, U.
%      SYS(5) = Reserved for root finding. Must be zero.
%      SYS(6) = Direct feedthrough flag (1=yes, 0=no). The s-function
%               has direct feedthrough if U is used during the FLAG=3
%               call. Setting this to 0 is akin to making a promise that
%               U will not be used during FLAG=3. If you break the promise
%               then unpredictable results will occur.
%      SYS(7) = Number of sample times. This is the number of rows in TS.
%
%
%      X0     = Initial state conditions or [] if no states.
%
%      STR    = State ordering strings which is generally specified as [].
%
%      TS     = An m-by-2 matrix containing the sample time
%               (period, offset) information. Where m = number of sample
%               times. The ordering of the sample times must be:
%
%               TS = [0      0,      : Continuous sample time.
%                     0      1,      : Continuous, but fixed in minor step
%                                      sample time.
%                     PERIOD OFFSET, : Discrete sample time where
%                                      PERIOD > 0 & OFFSET < PERIOD.
%                     -2     0];     : Variable step discrete sample time
%                                      where FLAG=4 is used to get time of
%                                      next hit.
%
%               There can be more than one sample time providing
%               they are ordered such that they are monotonically
%               increasing. Only the needed sample times should be
%               specified in TS. When specifying more than one
%               sample time, you must check for sample hits explicitly by
%               seeing if
%                  abs(round((T-OFFSET)/PERIOD) - (T-OFFSET)/PERIOD)
%               is within a specified tolerance, generally 1e-8. This
%               tolerance is dependent upon your model's sampling times
%               and simulation time.
%
%               You can also specify that the sample time of the S-function
%               is inherited from the driving block. For functions which
%               change during minor steps, this is done by
%               specifying SYS(7) = 1 and TS = [-1 0]. For functions which
%               are held during minor steps, this is done by specifying
%               SYS(7) = 1 and TS = [-1 1].
%
%      SIMSTATECOMPLIANCE = Specifices how to handle this block when saving and
%                           restoring the complete simulation state of the
%                           model. The allowed values are: 'DefaultSimState',
%                           'HasNoSimState' or 'DisallowSimState'. If this value
%                           is not speficified, then the block's compliance with
%                           simState feature is set to 'UknownSimState'.

switch flag,
  case 0,
    [sys,x0,str,ts]=mdlInitializeSizes;
  case 1,
    sys=mdlDerivatives(t,x,u);
  case {2,4,9},
    sys=[];
  case 3,
    sys=mdlOutputs(t,x,u);
  otherwise
    error(['Unhandled flag = ',num2str(flag)]);
end

%
%=============================================================================
% mdlInitializeSizes
% Return the sizes, initial conditions, and sample times for the S-function.
%=============================================================================
%
function [sys,x0,str,ts]=mdlInitializeSizes

%
% call simsizes for a sizes structure, fill it in and convert it to a
% sizes array.
%
% Note that in this example, the values are hard coded.  This is not a
% recommended practice as the characteristics of the block are typically
% defined by the S-function parameters.
%
sizes = simsizes;

sizes.NumContStates  = 9;
sizes.NumDiscStates  = 0;
sizes.NumOutputs     = 9;
sizes.NumInputs      = 9;
sizes.DirFeedthrough = 0;
sizes.NumSampleTimes = 1;   % at least one sample time is needed

sys = simsizes(sizes);

%
% initialize the initial conditions
%
x0  = [0 0 0 0 0 0 0 0 0];

%
% str is always an empty matrix
%
str = [];

%
% initialize the array of sample times
%
ts  = [0 0];
glo;
% end mdlInitializeSizes

%
%=============================================================================
% mdlDerivatives
% Return the derivatives for the continuous states.
%=============================================================================
%
function sys=mdlDerivatives(t,x,u)
global kt km 
global Lap m p
global vin Vb ns
global y
ne1=u(1);ne2=u(2);ne3=u(3);ne4=u(4);
Fi1=u(5);Fi2=u(6);Fi3=u(7);Fi4=u(8);
ns=u(9);

pb=x(1);qb=x(2);rb=x(3);
ue=x(4);ve=x(5);we=x(6);
fi=x(7);the=x(8);pus=x(9);

%转动惯量
Ixx=0.26;
Iyy=0.32;
Izz=0.35;

g=9.8;

%坐标转换
Tbe=[1 0 0;0 cos(fi) sin(fi);0 -sin(fi) cos(fi);]*[cos(the) 0 -sin(the);...
    0 1 0;sin(the) 0 cos(the);]*[cos(fi) sin(fi) 0;-sin(fi) cos(fi) 0;0 0 1;];
Trub=[cos(pi/4) -sin(pi/4) 0;sin(pi/4) cos(pi/4) 0;0 0 1];

Vb=Tbe*[ue ve we]';
ub=Vb(1);vb=Vb(2);wb=Vb(3);

%主桨输出
A=pi*0.35^2;
Tmp=3*quad(@mpropellert,0.08,0.32);
Mmp=3*quad(@mpropellerm,0.08,0.31);
vin=-wb/2+((wb/2)^2+Tmp/(2*p*A))^0.5;
%副桨输出
%副桨力矩计算
%暂用效率代替转速，取kt=4*9.8，km=0.2*9.8
Tap1=kt*ne1^2;
Tap2=kt*ne2^2;
Tap3=kt*ne3^2;
Tap4=kt*ne4^2;
Tap=kt*(ne1^2+ne2^2+ne3^2+ne4^2);
Map=[Lap*sin(pi/4)*(Tap2+Tap3-Tap1-Tap4) Lap*sin(pi/4)*(Tap3+Tap4-Tap1-Tap2)...
    km*(ne1^2+ne3^2-ne2^2-ne4^2)]';

%舵片输出
Omg=[pb,qb,rb]';
[Fru,Mru]=rudder(Trub*Omg,Trub*Vb,Fi1,Fi2,Fi3,Fi4);
Fru=Trub\Fru;
Mru=Trub\Mru;

%机身
%[Fbody,Mbody]=body(Vb);
S=0.5;
Fbody=-0.5*p*S*((wb+vin)^2+ub^2+vb^2)^0.5*[ub vb wb+vin]'*0.1;
%重力
Fg=Tbe*m*[0 0 -g]';
%合力计算
F=Fru+Fbody+Fg+[0 0 Tap]'+[0 0 Tmp]';
%M=Mru+Map+[0 0 Mmp]'; 
M=Map;

y(1)=M(1)/Ixx-qb*rb*(Izz-Iyy)/Ixx;
y(2)=M(2)/Iyy-pb*rb*(Ixx-Izz)/Iyy;
y(3)=M(3)/Izz-pb*qb*(Iyy-Ixx)/Izz;

invub=F(1)/m+vb*rb-wb*qb;
invvb=F(2)/m+wb*pb-ub*rb;
invwb=F(3)/m+ub*qb-vb*pb;

y(4:6)=Tbe\[invub invvb invwb]';

y(7:9)=[1,tan(the)*sin(fi),cos(fi)*tan(the);
    0,cos(fi),-sin(fi);
    0,sin(fi)/cos(the),cos(fi)/cos(the)]*[pb qb rb]';
sys = y;
% end mdlDerivatives


%
%=============================================================================
% mdlOutputs
% Return the block outputs.
%=============================================================================
%
function sys=mdlOutputs(t,x,u)
global y
sys = y;

% end mdlOutputs


function glo
global l L Lap vin p  kt km m
%结构参数
l=0.23;
L=0.12;
Lap=0.55;

%环境参数
p=1.292;
m=15;

%电机系数
kt=4*9.8;
km=0.2*9.8;

vin=10;


function dM=mpropellerm(r)
global ns Vb vin p
% # 初值预设
%角速度转换
b=0.04;
wb=Vb(3);
Rig=15-28*(r-0.08);
%Rig=acot(5.2*r-0.416)/7.2*180/pi;
% # 微元分析
W=((wb+vin)^2+(2*pi*ns*r).^2).^0.5;
bet=atan((wb+vin)./(2*pi*ns*r));
alf=Rig-bet*180/pi;

[cl,cd]=NACA0008(alf);

dD=0.5*b*p*cd.*W.^2;
dL=0.5*b*p*cl.*W.^2;

dM=r.*(dL.*sin(bet)+dD.*cos(bet));


function dT=mpropellert(r)
global ns Vb vin p
% # 初值预设
%角速度转换
b=0.04;
wb=Vb(3);
Rig=15-28*(r-0.08);
%Rig=acot(5.2*r-0.416)/7.2*180/pi;
%Rig=acot(6.4*r+0.488)/3;
% # 微元分析
W=((wb+vin)^2+(2*pi*ns*r).^2).^0.5;
bet=atan((wb+vin)./(2*pi*ns*r));
alf=Rig-bet*180/pi;

[cl,cd]=CLARK(alf);

dD=0.5*b*p*cd.*W.^2;
dL=0.5*b*p*cl.*W.^2;

dT=dL.*cos(bet)-dD.*sin(bet);


function [Fru,Mru]=rudder(Omg,Vb,Fi1,Fi2,Fi3,Fi4)
global p vin 
global l L
S=0.0265;

T1=[0 L 0; L 0 l; 0 -l 0];
T2=[0 L -l;L 0 0;l 0 0];
T3=[0 L 0;L 0 -l;0 l 0];
T4=[0 L l;L 0 0;-l 0 0];

V1=T1*Omg+Vb;
V2=T2*Omg+Vb;
V3=T3*Omg+Vb;
V4=T4*Omg+Vb;

u1=V1(1);v1=V1(2);w1=V1(3);
u2=V2(1);v2=V2(2);w2=V2(3);
u3=V3(1);v3=V3(2);w3=V3(3);
u4=V4(1);v4=V4(2);w4=V4(3);

bet1=atan(v1/(vin+w1));
bet3=-atan(v3/(vin+w3));
bet2=-atan(u2/(vin+w2));
bet4=atan(u4/(vin+w4));

alf1=Fi1+bet1*180/pi;
alf2=Fi2+bet2*180/pi;
alf3=Fi3+bet3*180/pi;
alf4=Fi4+bet4*180/pi;

[cl1,cd1]=CLARKy(alf1);
[cl2,cd2]=CLARKy(alf2);
[cl3,cd3]=CLARKy(alf3);
[cl4,cd4]=CLARKy(alf4);

V1abs=(v1^2+(w1+vin)^2)^0.5;
V3abs=(v3^2+(w3+vin)^2)^0.5;
V2abs=(u2^2+(w2+vin)^2)^0.5;
V4abs=(u4^2+(w4+vin)^2)^0.5;

D1=0.5*p*S*V1abs^2*cd1;
D2=0.5*p*S*V2abs^2*cd2;
D3=0.5*p*S*V3abs^2*cd3;
D4=0.5*p*S*V4abs^2*cd4;

L1=0.5*p*S*V1abs^2*cl1;
L2=0.5*p*S*V2abs^2*cl2;
L3=0.5*p*S*V3abs^2*cl3;
L4=0.5*p*S*V4abs^2*cl4;

F1=L1*cos(bet1)+D1*sin(bet1);
F2=L2*cos(bet2)+D2*sin(bet2);
F3=L3*cos(bet3)+D3*sin(bet3);
F4=L4*cos(bet4)+D4*sin(bet4);

Q1=L1*sin(bet1)-D1*cos(bet1);
Q2=L2*sin(bet2)-D2*cos(bet2);
Q3=L3*sin(bet3)-D3*cos(bet3);
Q4=L4*sin(bet4)-D4*cos(bet4);

Fru=[F2-F4 -F1+F3 Q1+Q2+Q3+Q4]';
Mru=[-L*(F1-F3)+l*(Q2-Q4) L*(-F2+F4)+l*(Q3-Q1) l*(-F2-F3-F1-F4)]';


function [cl,cd]=CLARKy(alf)
%ALF=-50:0.5:50;
ALF=[-50 -49.5 -49 -47 -46.5 -46 -45.5 -45 -44.5 -44 -43.5 -43 -42.5...
     -41 -40.5 -40 -39.5 -39 -38.5 -38 -37.5 -36 -35.5 -35 -34.5 -34 ...
     -33.5 -33 -32 -31.5 -31 -30.5 -30 -29.5 -28.5 -28 -27.5 -27 -26.5...
     -26 -25.5 -25 -24.5 -24 -23.5 -22.5 -22 -21.5 -21 -20 -19.5 -19 ...
     -18.5 -18 -17.5 -17 -16.5 -16 -15.5 -15 -14.5 -14 -13.5 -13 -12.5...
     -12 -11.5 -11 -10.5 -10 -9.5 -9 -8.5 -8 -7.5 -7 -6.5 -6 -5.5 -5 ...
     -4.5 -4 -3.5 -3 -2.5 -2 -1.5 -1 -0.5 0 0.5 1 1.5 2 2.5	3 3.5 4 4.5...
     5 5.5 6 6.5 7 7.5 8 8.5 9 9.5 10 10.5 11 11.5 12 12.5 13 13.5 14 ...
     14.5 15 15.5 16 16.5 17 17.5 18 18.5 19.5 20 20.5 21 21.5 22 22.5...
     23	23.5 24	24.5 25	25.5 26	26.5 27	27.5 28	28.5 29	29.5 30	30.5...
     31	31.5 32	32.5 33	33.5 34	34.5 35	35.5 36.5 37 37.5 38 38.5 39.5...
     40	40.5 41	41.5 42.5 43 43.5 44 44.5 45 46	46.5 47	47.5 48	48.5...
     49.5 50];
data_CD=[0.4658 0.4639 0.4618 0.4519 0.4490 0.4462 0.4432 0.4401...
         0.4368 0.4334 0.4298 0.4261 0.4222 0.4100 0.4057 0.4015...
         0.3971 0.3927 0.3881 0.3835 0.3787 0.3641 0.3592 0.3543...
         0.3493 0.3443 0.3392 0.3342 0.3236 0.3184 0.3133 0.3082...
         0.3030 0.2980 0.2874 0.2822 0.2771 0.2720 0.2671 0.2632...
         0.2576 0.2517 0.2467 0.2419 0.2375 0.2276 0.2224 0.2177...
         0.2141 0.2040 0.1993 0.1951 0.1958 0.1865 0.1818 0.1781...
         0.1807 0.1697 0.1657 0.1652 0.1689 0.1541 0.1503 0.1508...
         0.1392 0.1342 0.1324 0.1218 0.1170 0.1100 0.1045 0.1018...
         0.0940 0.0911 0.0894 0.0832 0.0807 0.0746 0.0682 0.0614...
         0.0358 0.0312 0.0279 0.0261 0.0245 0.0229 0.0208 0.0194...
         0.0190 0.0188 0.0185 0.0182 0.0181 0.0178 0.0176 0.0173...
         0.0172 0.0173 0.0175 0.0180 0.0186 0.0194 0.0201 0.0209...
         0.0219 0.0229 0.0240 0.0253 0.0267 0.0280 0.0295 0.0315...
         0.0339 0.0367 0.0401 0.0444 0.0507 0.0587 0.0662 0.0729...
         0.0799 0.0868 0.0954 0.1052 0.1176 0.1407 0.1820 0.2100...
         0.2151 0.2238 0.2273 0.2331 0.2433 0.2456 0.2513 0.2587...
         0.2648 0.2701 0.2759 0.2850 0.2888 0.2949 0.3009 0.3095...
         0.3139 0.3200 0.3262 0.3334 0.3392 0.3454 0.3513 0.3590...
         0.3643 0.3705 0.3764 0.3829 0.3891 0.3953 0.4014 0.4071...
         0.4190 0.4251 0.4308 0.4364 0.4416 0.4528 0.4584 0.4637...
         0.4686 0.4733 0.4831 0.4881 0.4927 0.4971 0.5011 0.5049...
         0.5130 0.5170 0.5206 0.5239 0.5269 0.5296 0.5356 0.5383...
        ];
data_CL=[-0.6143 -0.6172 -0.6198 -0.6281 -0.6299 -0.6313 -0.6325...
         -0.6333 -0.6340 -0.6344 -0.6345 -0.6344 -0.6340 -0.6315...
         -0.6303 -0.6287 -0.6269 -0.6248 -0.6224 -0.6198 -0.6168...
         -0.6068 -0.6029 -0.5987 -0.5943 -0.5896 -0.5846 -0.5795...
         -0.5687 -0.5627 -0.5565 -0.5500 -0.5434 -0.5368 -0.5229...
         -0.5151 -0.5073 -0.4995 -0.4916 -0.4852 -0.4772 -0.4669...
         -0.4579 -0.4492 -0.4413 -0.4245 -0.4137 -0.4042 -0.3974...
         -0.3789 -0.3682 -0.3597 -0.3649 -0.3453 -0.3341 -0.3274...
         -0.3414 -0.3147 -0.3059 -0.3128 -0.3429 -0.2929 -0.2920...
         -0.3252 -0.2845 -0.2819 -0.3212 -0.2755 -0.2857 -0.2830...
         -0.2774 -0.3384 -0.2974 -0.3402 -0.4137 -0.3818 -0.4669...
         -0.4851 -0.4808 -0.4215 -0.3330 -0.2533 -0.1808 -0.1075...
         -0.0348 0.0369 0.1062 0.2164 0.2951 0.3647 0.4308 0.4955...
         0.5478 0.6065 0.6610 0.7200 0.7725 0.8212 0.8729 0.9215...
         0.9703 1.0184 1.0644  1.1090 1.1536 1.1976 1.2396 1.2762...
         1.3041 1.3255 1.3366 1.3477 1.3596 1.3687 1.3697 1.3646...
         1.3458 1.3198 1.3051 1.3000 1.2953 1.2929 1.2830 1.2675...
         1.2413 1.1730 0.7454 0.7073 0.6961 0.7256 0.7031 0.7080...
         0.7363 0.7153 0.7202 0.7363 0.7339 0.7320 0.7392 0.7532...
         0.7456 0.7495 0.7564 0.7660 0.7613 0.7651 0.7716 0.7753...
         0.7748 0.7778 0.7817 0.7877 0.7856 0.7876 0.7900 0.7939...
         0.7935 0.7941 0.7951 0.7962 0.7972 0.7970 0.7970 0.7967...
         0.7963 0.7942 0.7930 0.7915 0.7898 0.7879 0.7831 0.7805...
         0.7778 0.7747 0.7714 0.7678 0.7600 0.7561 0.7519 0.7473...
         0.7426 0.7375 0.7271 0.7218];  
cd=interp1(ALF,data_CD,alf,'pchip');
cl=interp1(ALF,data_CL,alf,'pchip');


function [cl,cd]=CLARK(alf)
ALF=-15:0.5:15;
data_CD=[0.09078 0.08795 0.08326 0.07812 0.07362 0.06993 0.06597...
         0.06270 0.05874 0.05499 0.05152 0.04882 0.04595 0.04336...
         0.04076 0.03816 0.03589 0.03341 0.03129 0.02891 0.02727...
         0.02536 0.02347 0.02177 0.01998 0.01367 0.01361 0.01357...
         0.01397 0.01403 0.01445 0.01342 0.01571 0.01593 0.01626...
         0.01582 0.01621 0.01707 0.01789 0.01820 0.01997 0.02124...
         0.02240 0.02362 0.02493 0.02656 0.02692 0.02979 0.03330...
         0.03827 0.04089 0.04098 0.06801 0.07193 0.07669 0.07918...
         0.08435 0.08824 0.09291 0.09893 0.10377 ];
data_CL=[-0.720 -0.698 -0.674 -0.649 -0.623 -0.595 -0.566 -0.534...
         -0.502 -0.469 -0.436 -0.402 -0.368 -0.333 -0.298 -0.263...
         -0.228 -0.192 -0.155 -0.119 -0.082 -0.045 -0.008 -0.133...
         -0.073 -0.013 0.044 0.101 0.157 0.214 0.227 0.274 0.336...
         0.388 0.439 0.488 0.539 0.586 0.633 0.683 0.729 0.774...
         0.819 0.859 0.902 0.941 0.983 1.007 1.023 1.012 1.040...
         1.089 0.912 0.943 0.972 0.998 1.023 1.048 1.071 1.092...
         1.110];
cd=interp1(ALF,data_CD,alf,'pchip');
cl=interp1(ALF,data_CL,alf,'pchip');


function [cl,cd]=NACA0008(alf)
ALF=-15:0.5:15;
data_CD=[0.16577 0.15275 0.13913 0.12845 0.11901 0.11011 0.10164 0.09372 0.08747 0.07996 0.07344 0.06906...
         0.06282 0.05868 0.05361 0.04965 0.04575 0.04214 0.03920 0.03605 0.03326 0.03075 0.02835 0.02608...
         0.01551 0.01526 0.01441 0.01425 0.01418 0.01416 0.01525 0.01410 0.01412 0.01419 0.01434 0.01519...
         0.01544 0.02575 0.02802 0.03042 0.03283 0.03572 0.03887 0.04182 0.04542 0.04932 0.05328 0.05835...
         0.06249 0.06873 0.07311 0.07963 0.08715 0.09339 0.10131 0.10979 0.11869 0.12813 0.13881 0.15243 0.16545];
data_CL=[-0.746 -0.749 -0.749 -0.746 -0.739 -0.730 -0.717 -0.700 -0.682 -0.661 -0.639 -0.614...
         -0.588 -0.559 -0.528 -0.495 -0.462 -0.428 -0.393 -0.357 -0.322 -0.286 -0.249 -0.212...
         -0.320 -0.265 -0.208 -0.153 -0.097 -0.042 0.000 0.057 0.113 0.170 0.226 0.282 0.338...
         0.216 0.253 0.290 0.326 0.361 0.397 0.431 0.466 0.499 0.532 0.562 0.591 0.618 0.643...
         0.665 0.685 0.704 0.720 0.733 0.742 0.749 0.752 0.752 0.749];
cd=interp1(ALF,data_CD,alf,'pchip');
cl=interp1(ALF,data_CL,alf,'pchip');










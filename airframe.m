function [sys,x0,str,ts] = airframe(t,x,u,flag)
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

sizes.NumContStates  = 12;
sizes.NumDiscStates  = 0;
sizes.NumOutputs     = 12;
sizes.NumInputs      = 2;
sizes.DirFeedthrough = 0;
sizes.NumSampleTimes = 1;   % at least one sample time is needed

sys = simsizes(sizes);

%
% initialize the initial conditions
%
global Pp Qp Rp Fip Sip Pup Uc Vc Wc X0 Y0 Z0 
x0  = [Pp Qp Rp Uc Vc Wc Fip Sip Pup X0 Y0 3000-Z0];

%
% str is always an empty matrix
%
str = [];

%
% initialize the array of sample times
%
ts  = [0 0];

% end mdlInitializeSizes

%
%=============================================================================
% mdlDerivatives
% Return the derivatives for the continuous states.
%=============================================================================
%
function sys=mdlDerivatives(t,x,u)
global wx wy wz;
global Sp Sb rp b c e mp mb Omg
global Pp Qp Rp Fip Sip Pup Uc Vc Wc X0 Y0 Z0 vz mode
deta=u(1);
dets=u(2);

g=9.5;                          %重力加速度
p=1.292;                        %空气密度
u=8.8*pi/180;                   %安装角

CDdet=0.0957;                   %气动阻力系数
CLdet=0.235;                    %气动升力系数
Clbet=-0.0014;                  %滚转力矩阻尼系数
Clp=-0.133;                     %滚转力矩阻尼系数
Clr=0.01;                       %滚转力矩阻尼系数
Cldet=-0.0063;                  %滚转力矩阻尼系数
Cmq=-1.864;                     %俯仰力矩阻尼系数
Cmdet=0.294;                    %俯仰力矩阻尼系数
Cnbet=0.0005;                   %偏航力矩阻尼系数
Cnp=-0.013;                     %偏航力矩阻尼系数
Cnr=-0.035;                     %偏航力矩阻尼系数
Cndet=0.0155;                   %偏航力矩阻尼系数
CYdet=0.1368;                   %侧向力系数
CYr=-0.006;                     %侧向力系数
CYbet=-0.0095;                  %侧向力系数

Pp=x(1);         
Qp=x(2);
Rp=x(3);                     %伞体绕三轴角速度
Uc=x(4);
Vc=x(5);
Wc=x(6);                        %节点在地理坐标系下速度
Fip=x(7);
Sip=x(8);
Pup=x(9);                       %伞体相对地理坐标系的欧拉角
X0=x(10);
Y0=x(11);
Z0=3000-x(12);

%附加质量计算
AR=b/c;                         %展弦比
Ka=0.85;                        %三维效应修正因子
Kb=1.0;                         %三维效应修正因子
af11=p*Ka*pi*e^2*b/4;
af22=p*Kb*pi*e^2*c/4;
af33=p*(AR/(AR+1))*pi*c^2*b/4;
af44=0.055*p*(AR/(AR+1))*b^3*c^2;
af55=0.0308*p*(AR/(AR+1))*c^4*b;
af66=0.055*p*b^3*e^2;           %平直翼下附加质量

zpc=rp*sin(Omg/2)/(Omg/2);
zrc=zpc*af22/(af22+af44/rp^2);
zpr=zpc-zrc;                    %滚转中心，俯仰中心距
h1=(Omg/2)/4;                    

ma11=(1+8/3*h1^2)*af11;
ma22=(4^2*af22+af44)/zpc^2;
ma33=af33;
Ia11=(zpr^2/zpc^2)*rp^2*af22+(zrc^2/zpc^2)*af44;
Ia22=af55;
Ia33=(1+8*h1^2)*af66;           %圆弧翼下附加质量矩阵

If=[ma11,0,0;0,ma22,0;0,0,ma33];%附加质量矩阵
Im=[Ia11,0,0;0,Ia22,0;0,0,Ia33];%附加转动惯量矩阵
%附加质量计算结束

Tp=[cos(Sip)*cos(Pup),cos(Sip)*sin(Pup),-sin(Sip);
    sin(Fip)*sin(Sip)*cos(Pup)-cos(Fip)*sin(Pup),sin(Fip)*sin(Sip)*sin(Pup)+cos(Fip)*cos(Pup),sin(Fip)*cos(Sip);
    cos(Fip)*sin(Sip)*cos(Pup)+sin(Fip)*sin(Pup),cos(Fip)*sin(Sip)*sin(Pup)-sin(Fip)*cos(Pup),cos(Fip)*cos(Sip)];
%地理坐标系到伞体坐标系坐标变换矩阵

WB=mb*g*[-sin(Sip);sin(Fip)*cos(Sip);cos(Sip)*cos(Fip)];        %负载重力
WP=(mp+mb)*g*[-sin(Sip);sin(Fip)*cos(Sip);cos(Sip)*cos(Fip)];   %负载伞体合重

Ip=[79.3,0,0;0,12.3,0;0,0,91];                                  %伞体转动惯量矩阵

Spw=[0,-Rp,Qp;Rp,0,-Pp;-Qp,Pp,0];                               %伞体转速叉乘矩阵
Scp=[0,-rp*cos(u),0;rp*cos(u),0,-rp*sin(u);0,rp*sin(u),0];      %伞体与节点距离叉乘矩阵

Sudup=Tp*[Uc;Vc;Wc]+Spw*[rp*sin(u);0;rp*cos(u)];                %伞体在伞体坐标系下速度
Up=Sudup(1);
Vp=Sudup(2);
Wp=Sudup(3);

Suduw=[wx;wy;wz];                                               %风速
Sudupw=Tp*Suduw;
Upw=Sudupw(1);Vpw=Sudupw(2);Wpw=Sudupw(3);

Vcabs=((Uc-Suduw(1))^2+(Vc-Suduw(2))^2+(Wc-Suduw(3))^2)^0.5;             %节点相对于风的速度值
Vpabs=((Up-Upw)^2+(Vp-Vpw)^2+(Vp-Vpw)^2)^0.5;                            %伞体相对于风的速度值
ap=atan((Wp-Wpw)/(Up-Upw));                                              %伞体攻角
bp=asin((Vp-Vpw)/Vpabs);                                                 %伞体侧滑角

ap1=ap*180/pi;
DELTAS=[0,0.5,1];
ALPHAP=linspace(-10,80,91);
if dets==0
data_CD=[linspace(0.12,0.15,11) linspace(0.1563,0.2,9) linspace(0.2118,0.33,11) linspace(0.3465,0.66,20)  linspace(0.676,0.82,10) linspace(0.823,0.85,10) linspace(0.849,0.84,10) linspace(0.836,0.8,10)];
data_CL=[0.35 0.355 0.36 0.368 0.378 0.385 0.395 0.405 0.42 0.435 0.45 linspace(0.473,0.68,10) linspace(0.6729,0.63,7) linspace(0.65,0.73,8) 0.715 0.7 0.67 0.63 0.58 linspace(0.574,0.55,5) linspace(0.552,0.56,5) linspace(0.558,0.55,5) linspace(0.544,0.46,15) linspace(0.4475,0.21,20)];
data_Cmc=[linspace(-0.06,-0.11,6) linspace(-0.1127,-0.15,15) linspace(-0.1572,-0.51,50) linspace(-0.5045,-0.4,20)];
CDP=interp1(ALPHAP,data_CD,ap1)+CDdet*abs(deta);
CL=interp1(ALPHAP,data_CL,ap1)+CLdet*abs(deta);
Cmc=interp1(ALPHAP,data_Cmc,ap1);
else
data_CD=[linspace(0.12,0.15,11) linspace(0.1563,0.2,9) linspace(0.2118,0.33,11) linspace(0.3465,0.66,20)  linspace(0.676,0.82,10) linspace(0.823,0.85,10) linspace(0.849,0.84,10) linspace(0.836,0.8,10);
        linspace(0.12,0.18,11) linspace(0.1965,0.84,40) linspace(0.85,0.94,10) linspace(0.946,1,10) linspace(1,1,10) linspace(0.995,0.95,10);
        linspace(0.12,0.26,11) linspace(0.276,0.9,40) linspace(0.905,1,20) linspace(0.997,0.97,10) linspace(0.961,0.88,10)];
data_CL=[0.35 0.355 0.36 0.368 0.378 0.385 0.395 0.405 0.42 0.435 0.45 linspace(0.473,0.68,10) linspace(0.6729,0.63,7) linspace(0.65,0.73,8) 0.715 0.7 0.67 0.63 0.58 linspace(0.574,0.55,5) linspace(0.552,0.56,5) linspace(0.558,0.55,5) linspace(0.544,0.46,15) linspace(0.4475,0.21,20);
        linspace(0.35,0.58,11) linspace(0.62,0.98,10) linspace(0.956,0.86,5) linspace(0.872,0.92,5) linspace(0.902,0.74,10) linspace(0.734,0.68,10) linspace(0.676,0.66,5) linspace(0.6467,0.46,15) linspace(0.451,0.37,10) linspace(0.354,0.21,10);
        linspace(0.43,1,18)+0.2 linspace(0.9825,0.86,8)+0.2 linspace(0.872,0.92,5)+0.2 linspace(0.892,0.78,5)+0.2 linspace(0.769,0.56,20)+0.2 linspace(0.544,0,35)+0.2];
data_Cmc=[linspace(-0.06,-0.11,6) linspace(-0.1127,-0.15,15) linspace(-0.1572,-0.51,50) linspace(-0.5045,-0.4,20);
        linspace(-0.06,-0.45,41) linspace(-0.454,-0.57,30) linspace(-0.571,-0.58,10) linspace(-0.574,-0.52,10);
        linspace(-0.11,-0.53,51) linspace(-0.533,-0.56,10) linspace(-0.561,-0.57,10) linspace(-0.5665,-0.5,20)];
CDP=interp2(ALPHAP,DELTAS,data_CD,ap1,dets)+CDdet*abs(deta);
CL=interp2(ALPHAP,DELTAS,data_CL,ap1,dets)+CLdet*abs(deta);
Cmc=interp2(ALPHAP,DELTAS,data_Cmc,ap1,dets);
end
%对气动力系数和俯仰力矩系数的实验值进行二维插值

Cx=(-CDP*(Up-Upw)+CL*(Wp-Wpw))/Vpabs;
Cy=CYbet*bp+CYr*Rp*b/2/Vpabs+CYdet*deta;
Cz=(-CDP*(Wp-Wpw)-CL*(Up-Upw))/Vpabs;                          %气动力在伞体坐标系下的系数

Cl=Clbet*bp+Clp*Pp*b/2/Vpabs+Clr*Rp*b/2/Vpabs+Cldet*deta;
Cm=Cmc+Cmq*Qp*c/2/Vpabs+Cmdet*abs(deta);
Cn=Cnbet*bp+Cnp*Pp*b/2/Vpabs+Cnr*Rp*b/2/Vpabs+Cndet*deta;%气动力矩在伞体坐标系下的系数

Fpa=0.5*p*(Vpabs)^2*Sp*[Cx;Cy;Cz];                       %气动力
Mpa=0.5*p*(Vpabs)^2*Sp*[b*Cl;c*Cm-0.25*c*Cz;b*Cn];       %气动力矩

Fp=WP+Fpa-mp*Spw*Spw*[rp*sin(u);0;rp*cos(u)]-Spw*If*Spw*[rp*sin(u);0;rp*cos(u)]-0.5*p*Vcabs*Tp*Sb*[Uc-Suduw(1);Vc-Suduw(2);Wc-Suduw(3)];%伞体受力
Mp=Mpa-Spw*(Ip+Im)*[Pp;Qp;Rp]-Scp*WB+Scp*0.5*p*Vcabs*Tp*Sb*[Uc-Suduw(1);Vc-Suduw(2);Wc-Suduw(3)];                                       %伞体所受力矩

B1=Mp;
B2=Fp;

Dx=[Ip+Im,-Scp*Tp*mb;
    -If*Scp-mp*Scp,If*Tp+(mp+mb)*Tp]\[B1;B2];                   %动力学方程
y(1:6)=Dx;

y(7:9)=[1,tan(Sip)*sin(Fip),cos(Fip)*tan(Sip);
    0,cos(Fip),-sin(Fip);
    0,sin(Fip)/cos(Sip),cos(Fip)/cos(Sip)]*[Pp;Qp;Rp];          %欧拉角求解
y(10:12)=[Uc;Vc;Wc];                                            %地理坐标系下位移

if mode==2
    y(12)=vz;
end

sys = y;
% end mdlDerivatives


%
%=============================================================================
% mdlOutputs
% Return the block outputs.
%=============================================================================
%
function sys=mdlOutputs(t,x,u)

sys = x;

% end mdlOutputs


















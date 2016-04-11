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

g=9.5;                          %�������ٶ�
p=1.292;                        %�����ܶ�
u=8.8*pi/180;                   %��װ��

CDdet=0.0957;                   %��������ϵ��
CLdet=0.235;                    %��������ϵ��
Clbet=-0.0014;                  %��ת��������ϵ��
Clp=-0.133;                     %��ת��������ϵ��
Clr=0.01;                       %��ת��������ϵ��
Cldet=-0.0063;                  %��ת��������ϵ��
Cmq=-1.864;                     %������������ϵ��
Cmdet=0.294;                    %������������ϵ��
Cnbet=0.0005;                   %ƫ����������ϵ��
Cnp=-0.013;                     %ƫ����������ϵ��
Cnr=-0.035;                     %ƫ����������ϵ��
Cndet=0.0155;                   %ƫ����������ϵ��
CYdet=0.1368;                   %������ϵ��
CYr=-0.006;                     %������ϵ��
CYbet=-0.0095;                  %������ϵ��

Pp=x(1);         
Qp=x(2);
Rp=x(3);                     %ɡ����������ٶ�
Uc=x(4);
Vc=x(5);
Wc=x(6);                        %�ڵ��ڵ�������ϵ���ٶ�
Fip=x(7);
Sip=x(8);
Pup=x(9);                       %ɡ����Ե�������ϵ��ŷ����
X0=x(10);
Y0=x(11);
Z0=3000-x(12);

%������������
AR=b/c;                         %չ�ұ�
Ka=0.85;                        %��άЧӦ��������
Kb=1.0;                         %��άЧӦ��������
af11=p*Ka*pi*e^2*b/4;
af22=p*Kb*pi*e^2*c/4;
af33=p*(AR/(AR+1))*pi*c^2*b/4;
af44=0.055*p*(AR/(AR+1))*b^3*c^2;
af55=0.0308*p*(AR/(AR+1))*c^4*b;
af66=0.055*p*b^3*e^2;           %ƽֱ���¸�������

zpc=rp*sin(Omg/2)/(Omg/2);
zrc=zpc*af22/(af22+af44/rp^2);
zpr=zpc-zrc;                    %��ת���ģ��������ľ�
h1=(Omg/2)/4;                    

ma11=(1+8/3*h1^2)*af11;
ma22=(4^2*af22+af44)/zpc^2;
ma33=af33;
Ia11=(zpr^2/zpc^2)*rp^2*af22+(zrc^2/zpc^2)*af44;
Ia22=af55;
Ia33=(1+8*h1^2)*af66;           %Բ�����¸�����������

If=[ma11,0,0;0,ma22,0;0,0,ma33];%������������
Im=[Ia11,0,0;0,Ia22,0;0,0,Ia33];%����ת����������
%���������������

Tp=[cos(Sip)*cos(Pup),cos(Sip)*sin(Pup),-sin(Sip);
    sin(Fip)*sin(Sip)*cos(Pup)-cos(Fip)*sin(Pup),sin(Fip)*sin(Sip)*sin(Pup)+cos(Fip)*cos(Pup),sin(Fip)*cos(Sip);
    cos(Fip)*sin(Sip)*cos(Pup)+sin(Fip)*sin(Pup),cos(Fip)*sin(Sip)*sin(Pup)-sin(Fip)*cos(Pup),cos(Fip)*cos(Sip)];
%��������ϵ��ɡ������ϵ����任����

WB=mb*g*[-sin(Sip);sin(Fip)*cos(Sip);cos(Sip)*cos(Fip)];        %��������
WP=(mp+mb)*g*[-sin(Sip);sin(Fip)*cos(Sip);cos(Sip)*cos(Fip)];   %����ɡ�����

Ip=[79.3,0,0;0,12.3,0;0,0,91];                                  %ɡ��ת����������

Spw=[0,-Rp,Qp;Rp,0,-Pp;-Qp,Pp,0];                               %ɡ��ת�ٲ�˾���
Scp=[0,-rp*cos(u),0;rp*cos(u),0,-rp*sin(u);0,rp*sin(u),0];      %ɡ����ڵ�����˾���

Sudup=Tp*[Uc;Vc;Wc]+Spw*[rp*sin(u);0;rp*cos(u)];                %ɡ����ɡ������ϵ���ٶ�
Up=Sudup(1);
Vp=Sudup(2);
Wp=Sudup(3);

Suduw=[wx;wy;wz];                                               %����
Sudupw=Tp*Suduw;
Upw=Sudupw(1);Vpw=Sudupw(2);Wpw=Sudupw(3);

Vcabs=((Uc-Suduw(1))^2+(Vc-Suduw(2))^2+(Wc-Suduw(3))^2)^0.5;             %�ڵ�����ڷ���ٶ�ֵ
Vpabs=((Up-Upw)^2+(Vp-Vpw)^2+(Vp-Vpw)^2)^0.5;                            %ɡ������ڷ���ٶ�ֵ
ap=atan((Wp-Wpw)/(Up-Upw));                                              %ɡ�幥��
bp=asin((Vp-Vpw)/Vpabs);                                                 %ɡ��໬��

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
%��������ϵ���͸�������ϵ����ʵ��ֵ���ж�ά��ֵ

Cx=(-CDP*(Up-Upw)+CL*(Wp-Wpw))/Vpabs;
Cy=CYbet*bp+CYr*Rp*b/2/Vpabs+CYdet*deta;
Cz=(-CDP*(Wp-Wpw)-CL*(Up-Upw))/Vpabs;                          %��������ɡ������ϵ�µ�ϵ��

Cl=Clbet*bp+Clp*Pp*b/2/Vpabs+Clr*Rp*b/2/Vpabs+Cldet*deta;
Cm=Cmc+Cmq*Qp*c/2/Vpabs+Cmdet*abs(deta);
Cn=Cnbet*bp+Cnp*Pp*b/2/Vpabs+Cnr*Rp*b/2/Vpabs+Cndet*deta;%����������ɡ������ϵ�µ�ϵ��

Fpa=0.5*p*(Vpabs)^2*Sp*[Cx;Cy;Cz];                       %������
Mpa=0.5*p*(Vpabs)^2*Sp*[b*Cl;c*Cm-0.25*c*Cz;b*Cn];       %��������

Fp=WP+Fpa-mp*Spw*Spw*[rp*sin(u);0;rp*cos(u)]-Spw*If*Spw*[rp*sin(u);0;rp*cos(u)]-0.5*p*Vcabs*Tp*Sb*[Uc-Suduw(1);Vc-Suduw(2);Wc-Suduw(3)];%ɡ������
Mp=Mpa-Spw*(Ip+Im)*[Pp;Qp;Rp]-Scp*WB+Scp*0.5*p*Vcabs*Tp*Sb*[Uc-Suduw(1);Vc-Suduw(2);Wc-Suduw(3)];                                       %ɡ����������

B1=Mp;
B2=Fp;

Dx=[Ip+Im,-Scp*Tp*mb;
    -If*Scp-mp*Scp,If*Tp+(mp+mb)*Tp]\[B1;B2];                   %����ѧ����
y(1:6)=Dx;

y(7:9)=[1,tan(Sip)*sin(Fip),cos(Fip)*tan(Sip);
    0,cos(Fip),-sin(Fip);
    0,sin(Fip)/cos(Sip),cos(Fip)/cos(Sip)]*[Pp;Qp;Rp];          %ŷ�������
y(10:12)=[Uc;Vc;Wc];                                            %��������ϵ��λ��

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


















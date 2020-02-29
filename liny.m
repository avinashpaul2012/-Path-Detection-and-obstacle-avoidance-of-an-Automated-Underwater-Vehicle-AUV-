
m = 2234.5;zb=-0.041;Iy=1700;Zw_dot=-500;W=21920.4;B=21898;Zww=-400;
Zqq=500;Mq_dot=-1692.3;
U=1;Zw_dot=-500;
Mw_dot=-2090.4;
Zuuds=500;
Muuds=500;

%linearization
x1=0;
syms  x2 x3 x4;
F=[%inv(m - Zw_dot)*((W-B)*cos(x4)+ Zww*x1*U + Zqq*x2*U + (m)*U*x2 + Zw_dot*x1);...
    inv(Iy-Mq_dot)*(zb*B*sin(x4) + Mq_dot*U*x2 + Mw_dot*x1 );...
    x1*cos(x4)-sin(x4)*U;...
    x2];
X=[x2 ;x3 ;x4];
G=[Zuuds*U^2;Muuds*U^2;0;0];

A=jacobian(F,X);
C=double(subs(A,[x2 x3 x4],[0 0 0]));

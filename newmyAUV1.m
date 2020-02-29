
function xdot =newmyAUV1(t,V)
u=V(1);
v=V(2);
w=V(3);
q=V(4);
r=V(5);
x=V(6);
y=V(7);
z=V(8);
phi=V(9);
theta=V(10);
sai=V(11);
p=0;
dr=0.2;
ds=0;

%AUV parameters.
m = 2234.5;zb=-0.041;Ix=700;Iy=1700;Iz=2000;
Xu_dot=-141.9;B=21898;W=21920.4;                     %W=m*g.
Xuu=-35.4;Yv_dot=-1715.4;Yr_dot=186.9;Yvv=-667.5;Yrr=-32.5;Nr_dot=-1349;Nrr=-310;
Nvv=433.8;Xvr=1715.4;Yur=103.4;Nur=-1427;Yuv=-346.76;Nuv=-686.08;Nv_dot=957;
Xwq =137;Xqq =9587.4;Xrr = 832;zq_dot=-1701.9;Mq_dot=-1692.3;Mw_dot=-2090.4;Xprop=80;
Yuudr=1000;
Zw_dot=200;
Zww=-400;
Zqq=500;
Zuuds=500;
Muuds=500;
Nuudr=1000;
%function begining

udot=inv(m-Xu_dot)*(-(W-B)*sin(theta) + Xuu*u*abs(u) + (Xwq - m)*w*q  + (Xqq)*q^2 + (Xvr + m)*v*r...
    + (Xrr)*r^2 + Xprop);

vdot=inv(m - Yv_dot )*((W-B)*cos(theta)*sin(phi) + Yvv*v*abs(v) + Yrr*r*abs(r) + (Yur - m)*u*r...
     +Yuv*u*v + Yuudr*u^2*dr);

wdot=inv(m - Zw_dot)*((W-B)*cos(theta)*cos(phi) + Zww*w*abs(w) + Zqq*q*abs(q) + (m)*u*q...
    + Zuuds*u^2*ds);

qdot=inv(Iy-Mq_dot)*(zb*B*sin(theta) +  ...
      Muuds*u^2*ds );

rdot =inv(Iz-Nr_dot)*(Nvv*v*abs(v) + Nrr*r*abs(r) + (Nur)*u*r....
    + Nuv*u*v + Nuudr*u^2*dr); %NHS=0

X1dot=[udot;vdot;wdot;qdot;rdot];

n1=[x;y;z];
J1n2=[cos(sai)*cos(theta) -sin(sai)*cos(phi)+cos(sai)*sin(theta)*sin(phi) sin(sai)*sin(theta)+cos(sai)*sin(theta)*cos(phi);...
    sin(sai)*cos(theta) cos(sai)*cos(phi)+sin(sai)*sin(theta)*sin(phi) -cos(sai)*sin(phi)+sin(sai)*sin(theta)*cos(phi);...
    -sin(theta) cos(theta)*sin(phi) cos(theta)*cos(phi)];
v1=[u;v;w];

n1dot=J1n2*v1;

n2=[phi;theta;sai];
J2n2=[1 sin(phi)*tan(theta) cos(phi)*tan(theta);0 cos(phi) -sin(phi);0 sin(phi)*inv(cos(theta)) cos(phi)*inv(cos(theta))];
v2=[p;q;r];

n2dot=J2n2*v2;

X2dot=[n1dot;n2dot];

xdot=[X1dot;X2dot];









 syms s
num=sym2poly(10^10*(5*s^2+7.4945*s+1.3235));
den=sym2poly(10^8*s^3+82800000*s^2+104508799*s+8711277); 
G=tf(num,den);
H=1;
kp=0.5;
kd=0;
ki=0;
C=pid(kp,ki,kd);
T=feedback(C*G,H);


tspan=0:0.01:800;
X0=[0 0 0 0 0 0 0 0 0 0 0];
[t,X]=ode45('newmyAUV1',tspan,X0);
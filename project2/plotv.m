% octave script file
x=0:0.001:1;
rho=10*sin(10*sin(pi*x)).*(cos(pi*x)).*(cos(pi*x))+sin(pi*x).*cos(10*sin(pi*x));
plot(x,rho);
figure(2);
v=sin(10*sin(pi*x))/(10*pi*pi)+x;
h=1./200.;
x_sd=h:h:199.*h;
v_sd=load('v1.txt');
plot(x,v,x_sd,v_sd,'o');


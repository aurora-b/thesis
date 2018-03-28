clear
%note here, we must use rk4 instead of ode45 because ode45 uses variable
%time step size.
tend=5;
g=9;
n=2^g; %grid points
b=2*pi; %length of x axis
delx= b/n; %width of space step
delt=0.1*delx;
w=length(0:delt:tend);
visc=delx^1.2
x= 0:delx:b-delx; %adds delx each time and specifies grid points
uinit=zeros(1,n); %preallocating u
for i=1:n
    uinit(i)= sin(x(i));
end
[t,u]=rk4(@ode45try2,0, delt, uinit, 1); %1 time step from 0 to delt


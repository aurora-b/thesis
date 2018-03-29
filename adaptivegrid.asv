clear
%note here, we must use rk4 instead of ode45 because ode45 uses variable
%time step size.
tend=5;
g=9;
n=2^g; %grid points
b=2*pi; %length of x axis
delx= b/n; %width of space step
delt=0.5*delx;
w=length(0:delt:tend);
visc=delx^1.2
x= 0:delx:b-delx; %adds delx each time and specifies grid points
uinit=zeros(1,n); %preallocating u
for i=1:n
    uinit(i)= sin(x(i));
end



len = length(uinit); %number of grid points
j=5; %make sure to change this in the file too. this should be such that 2^(j-1) = len
lev   = 8;
nbcol = 100;
App=zeros(lev, len/2); 
Dt=zeros(lev,len/2);  
eps=5E-3;

%perform decomposition
[App(1,1:len/2),Dt(1,1:len/2)]=waveinter(uinit,1,0);
for i=2:lev
     Ex = App(i-1,1:(len)/(2^(i-1)));
    [App(i,1:((len/(2^i)))),Dt(i,1:(len/(2^i)))] = waveinter(Ex, 1,0);
end
[App,Dt,y1]=activegrid(App,Dt, uinit, eps, lev);

%% Find the active grid at EACH level.
%Restructure data structures so they are all 'len' in length
Dt1=zeros(lev,len);
App1=zeros(lev,len);
acgrid=zeros(lev+1,len);
for i=1:lev
    x=Dt(i,1:len/(2^i));
    Dt1(i,(2^i):(2^i):end)=x;
    x1=App(i,1:len/(2^i));
    App1(i,2^(i-1):(2^i):end)=x1;
    acgrid(i+1,:)=Dt1(i,:)+App1(i,:);
end
acgrid(1,:)=y1
%we now have an active grid structure. the first level is the finest level,
%and the next levels down go coarser and coarser

%%







[t,u]=rk4(@ode45try2,0, delt, uinit, 1); %1 time step from 0 to delt

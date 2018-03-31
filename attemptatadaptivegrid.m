%% ALL THIS CODE FOLLOWS FROM THE WORK DONE ON PURPLE PAPER ON 22 MARCH
%This doesn't implement solving on the adaptive grid, but shows where the
%adaptive grid should be based on active wavelets/perfect
%reconstruction/security zone

%% To setup this function, either run the odetry3 function, or run the below chunk.
%This chunk will reference ode45try2 and the rk4 function
clear
tend=6.1;
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
[t,u]=rk4(@ode45try2,0, tend, uinit, w); 
%%
for n=1:w
s = [u(n,:)]; 
len = length(s); %number of grid points
lev   = 8;
nbcol = 100;
App=zeros(lev, (length(s))/2); 
Dt=zeros(lev,[length(s)]/2);  
eps=5E-3;

%perform decomposition
[App(1,1:len/2),Dt(1,1:len/2)]=waveinter(s,1,0);
for i=2:lev
     Ex = App(i-1,1:(len)/(2^(i-1)));
    [App(i,1:((len/(2^i)))),Dt(i,1:(len/(2^i)))] = waveinter(Ex, 1,0);
end
[App, Dt]=activegrid(App,Dt, u(1,:), eps, lev);

%NOW we have data structures containing all the significant grid points
%based on the significant wavelet coefficients, the perfect reconstruction
%zone, and the coarsest level of scaling function decomposition

%Restructure data structures so they are all 'len' in length

%Restructure data structures so they are all 'len' in length
Dt1=zeros(lev,len);
App1=zeros(lev,len);
acgrid=zeros(lev,len);
for i=1:lev
    x=Dt(i,1:len/(2^i));
    Dt1(i,(2^i):(2^i):end)=x;
    x1=App(i,1:len/(2^i));
    App1(i,2^(i-1):(2^i):end)=x1;
    acgrid(i+1,:)=Dt1(i,:)+App1(i,:);
end
%acgrid(1,:)=y1;
%we now have an active grid structure. the first level is the finest level,
%and the next levels down go coarser and coarser

%now we sum all of these together to see where we have points that aren't 0
%(i.e. where we should preserve the grid points)
agrid=sum(acgrid,1);

% This chunk plots the active grid on the x axis
J=find(abs(agrid==0));
agrid(J)=NaN(size(J)); %set things below threshold to NaN so they don't plot

J1=find(abs(agrid)>=0);
agrid(J1)=-ones(size(J1));

%END of the chunk which plots the active grid on the x axis
figure(1)
plot(s);
hold on;
plot(agrid,'.','MarkerSize', 10);
axis([0 len -1 1])
set(gca, 'XTick', [0:0.1:1]*len, 'XTickLabel', [0:0.1:1]*2)
hold off;
mov(n)=getframe(figure(1));
end
vv = VideoWriter('activegrid_delt0.1timesdelx_gridpoints2power9_viscdelxpowerof1.2run2.avi');
vv.FrameRate = 110;  % Default 30
vv.Quality = 100;    % Default 75
open(vv)
writeVideo(vv,mov)
close(vv)  



  
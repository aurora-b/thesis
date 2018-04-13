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
[App,Dt]=activegrid(App,Dt, eps, lev);

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
%acgrid(1,:)=y1
%we now have an active grid structure. the first level is the finest level,
%and the next levels down go coarser and coarser

%% applying the numerical solution on the adaptive grid
% for i=1:lev
%    if i=1
%        yt=waveinterinv(App(1,:),Dt(1,:),1)
%        I3=0; %initialize a vector
%        I1=find(abs(Dt(1,:)>0);
%        I2=find(abs(App(1,:)>0); %find where the points are that weren't set to 0
%        for k=1:len/2
%        if ismember(k,I1)==1 %adding in neighbours of any significant wavelet point
%            I2(length(I2)+1)=k;
%            if k=len/2
%                 I2(length(I2)+1)=1
%            else
%                I2(length(I2)+1)=k+1
%            end
%        end
%        if ismember(k,I2)==1 %adding in neighbours of the scaling function points that we want to keep at each level
%            I1(length(I1)+1)=k;
%            if k=1
%                I1(length(I2)+1)=len/2
%            else
%                I1(length(I2)+1)=k-1
%            end
%        end
%        %now, we need to put this as a vector combined (find where on yt, we
%        %want to keep pts)
%        for ismember(k,I2)==1
%            I3(length(I3)+1)=2*k-1;
%        end
%        for ismember(k,I1)==1
%            I3(length(I3)+1)=2*k;
%        end
%        end 
%        %find consecutive chunks
%        I3=I3(2:end); %cut out the initialize value 0
%        I4=find(diff(I3)~=1) %if a value in this vector after 1 is n, then take 1:n as a vector
%        
%        
%        
%    end
%     
%     
%     
%     
% end
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% %%
% [t,u]=rk4(@ode45try2,0, delt, uinit, 1); %1 time step from 0 to delt

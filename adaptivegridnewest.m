%% WE ADD IN GRID POINTS ABOVE ACTIVE WAVELETS. SHOULD WE?
%% This code provides a way to create an adaptive grid! yay. this code depends on activegridcalc.m (which depends on activegrid.m), rk4setup.m, rk4try2.m, waveinter.m and waveinterinv.m
clear
tend=5;
g=9;
n=2^g; %grid points
b=2*pi; %length of x axis
delx= b/n; %width of space step
delt=10*delx;
w=length(0:delt:tend);
visc=delx^2/8;
x= 0:delx:b-delx; %adds delx each time and specifies grid points
uinit=zeros(1,n); %preallocating u
for i=1:n
    uinit(i)= sin(x(i));
end
%%
s = [uinit]; 
len = length(s); %number of grid points
lev   = 6;
yfd=zeros(lev,len,w+1);
for p=1:w
App=zeros(lev, (length(s))/2); 
Dt=zeros(lev,[length(s)]/2);  
eps=1E-5;
%perform decomposition
[App(1,1:len/2),Dt(1,1:len/2)]=waveinter(s,1,0);
for i=2:lev
     Ex = App(i-1,1:(len)/(2^(i-1)));
    [App(i,1:((len/(2^i)))),Dt(i,1:(len/(2^i)))] = waveinter(Ex, 1,0);
end
[App, Dt,tr]=activegridcalc(App, Dt, eps, lev); %get the app and dt matrices that we can perform finite difference on. also get the tracker matrix
%undo lifting
yt= zeros(lev, (length(x)));
yt(1, 1:((len/(2^(lev-1)))))=waveinterinv(App(lev,1:((len/(2^lev)))),Dt(lev,1:(len/(2^lev))),1); %start with y1, reconstruct
for i=2:lev
    yt(i, 1:(len/(2^(lev-i))))=waveinterinv(yt(i-1, 1:(len/(2^((lev-i+1))))),Dt((lev+1-i),1:(len/(2^(lev-i+1)))),1); %reconstruct up all levels
end
yt=flipud(yt);
for r=1:lev
yfd(r,1:(2^(r-1)):end,1)=yt(r,1:len/(2^(r-1)));
end

[t,thing]=rk4try2(@rk4setup,delt*(p-1), delt*p, yfd(1,:,p), 1,len,1); %go up one time step only and on the highest level
yfd(1,:,p+1)=thing(2,:); %found the next time step at the finest level
yfd(1,:,p+1)=yfd(1,:,p+1).*tr(1,:); %zeroes out the non significant points



%now find if we can restrict, if not then do finite difference/volume 
for r=2:lev
QQ=zeros(1,len);
for k=1:len
    Q=0;
    if tr(r,k)==1
        for f=1:r-1
            Q(length(Q)+1)=tr(f,k); %store the stuff above
        end 
    if ismember(1,Q)==1 %if there is a 1 in there
   QQ(k)= min(find(Q==1)-1);%if one of the values above ours is a 1, store in the kth spot of the QQ matrix the finest row that we have a value at
    else 
   QQ(k)=0; %if there are no singificant ones above it, print a 0 in the kth spot
    end
    end
end

if length(find(QQ~=0))==length(find(tr(r,:)==1)) %if every significant value on row r has one above it
I3=find(QQ~=0); %find indices where there are levels above with significant pts
for k=1:len
if ismember(k,I3)==1
    yfd(r,k,p)=yfd(QQ(k),k,p); %replace those points with the solution at the finer level
end
end
yfd(r,:,p)=yfd(r,:,p).*tr(r,:); %zeroes out the non significant points

else
%calculate finite-volume/finite-difference scheme on row r
[t,thing]=rk4try2(@rk4setup,delt*(p-1), delt*p, yfd(r,1:(2^(r-1)):end,p), 1,len,r); %go up one time step only and on the highest level
yfd(r,1:(2^(r-1)):end,p+1)=thing(2,:); %found the next time step at the finest level
yfd(r,:,p+1)=yfd(r,:,p+1).*tr(r,:); %zeroes out the non significant points

%but if there are any values on finer levels, we need to replace them
%because that will make this more accurate
I3=find(QQ~=0); %find indices where there are levels above with significant pts
for k=1:len
if ismember(k,I3)==1
    yfd(r,k,p)=yfd(QQ(k),k,p); %replace those points with the solution at the finer level
end
end
end
end
end

y=yfd(1,:,p); %start with finest level
%if any other level in a certain place isn't NaN, fill in that value
for k=1:len
for r=2:lev
if isnan(y(k))
    y(k)=yfd(r,k,p);
end  
end
end
s=y;



%plot it even though it has NaN in it
I = ~isnan(y);
plot(x(I),y(I))

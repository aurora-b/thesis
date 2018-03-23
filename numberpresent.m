for n=1:j
s = [u(n,:)]; 
len = length(s); %number of grid points

lev   = 4;
nbcol = 100;
App=zeros(lev, (length(s))/2); 
Dt=zeros(lev,[length(s)]/2); 
eps=1E-5

%perform decomposition
[App(1,:),Dt(1,:)]=waveinter(s,1,0);

for i=2:lev
     Ex = App(i-1,1:((length(x))/(2^(i-1))));
    [App(i,1:((len/(2^i)))),Dt(i,1:(len/(2^i)))] = waveinter(Ex, 1,0);
end
I2=find(abs(Dt)>eps);
numpres(n)=prod(size(I2));



end
plot(numpres) %what we notice is for large eps, numpres increases over time. For small eps, numpres decreases over time
xlabel('Time step')
ylabel('Number of Collocation Points above Threshold')
  
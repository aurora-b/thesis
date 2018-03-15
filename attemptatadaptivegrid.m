for n=1:j
s = [u(n,:)]; 
len = length(s); %number of grid points

lev   = 4;
nbcol = 100;
App=zeros(lev, (length(s))/2); 
Dt=zeros(lev,[length(s)]/2); 
cfd=zeros(lev,[length(s)]);

%perform decomposition
[App(1,:),Dt(1,:)]=waveinter(s,1,0);

for i=2:lev
     Ex = App(i-1,1:((length(x))/(2^(i-1))));
    [App(i,1:((len/(2^i)))),Dt(i,1:(len/(2^i)))] = waveinter(Ex, 1,0);
end

end

  
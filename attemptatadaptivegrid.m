%% ALL THIS CODE FOLLOWS FROM THE WORK DONE ON PURPLE PAPER ON 22 MARCH
%We start working at the first time step. This is step 0 in Rong
r=length(u(:,1))
n=1 
s = [u(n,:)]; 
len = length(s); %number of grid points
j=5; %make sure to change this in the file too. this should be such that 2^(j-1) = len

lev   = 3;
nbcol = 100;
App=zeros(lev, (length(s))/2); 
Dt=zeros(lev,[length(s)]/2);  
eps=1E-5

%perform decomposition
[App(1,1:len/2),Dt(1,1:len/2)]=waveinter(s,1,0);

for i=2:lev
     Ex = App(i-1,1:((length(x))/(2^(i-1))));
    [App(i,1:((len/(2^i)))),Dt(i,1:(len/(2^i)))] = waveinter(Ex, 1,0);
end
I2=find(abs(Dt')>eps); %finds points on the grid which are currently above the threshold. counts down column usually, but we want it to count down row.

%Now we need to find grid points which are adjacent in both scale and
%space. This is step 1 in Rong

%ADJACENT IN SPACE

%define which pts are the first pt on each scale
for i=1:lev
    firstpt(i)=(i-1)*2^(j-1)+1;
end

%define which pts are the last pt on each scale
for i=1:lev
    lastpt(i)=(i-1)*2^(j-1)+2^(j-i);
end

I3=I2;%store values in I3
for k=1:len*lev
     for i=1:lev
        if k>(len/2)*(i-1) && k<=(len/2)*(i) %determine which line the point is on.
            row=i;
        end
    end
    if ismember(k,I2)==1 && ismember(k,firstpt)==0 && ismember(k,lastpt)==0 %if k is significant but is not a first or last pt at each scale
        I3(length(I3)+1)=k-1; %add into the I2 array the wavelet prior in space
        I3(length(I3)+1)=k+1; %add into the I2 array the wavelet after in space
    end
    if ismember(k,I2)==1 && ismember(k,firstpt)==1 %if is the first point on some scale
        I3(length(I3)+1)=k+1; %add the second point on that scale
        I3(length(I3)+1)=k+2^(j-row)-1; %add the last point on that scale
    end
    if ismember(k,I2)==1 && ismember(k,lastpt)==1
        I3(length(I3)+1)=k-1;
        I3(length(I3)+1)=k-2^(j-row)+1;
    end
end


%ADJACENT IN SCALEI
%We want to add indices to the I3 matrix which correlate to be adjacent in
%scale.
for k=1:len*lev
if ismember(k,I2)==1 %if we are in the I2 matrix
    for i=1:lev
        if k>(len/2)*(i-1) && k<=(len/2)*(i) %determine which line the point is on.
            row=i;
        end
    end
        m=k-(len/2)*(row-1) %gives position on line
        if mod(m,2)==0 && row~=1 && row~=lev %if an even point on the line, and the line isn't the first or last
            I3(length(I3)+1)=k+(len/2-m) + m/2; %goes down one line
            I3(length(I3)+1)=k-m-(len/2)+2*m; % goes up one line
        end
        if mod(m,2)~=0 && row~=1 && row~=lev || row==lev %if odd on a middle line, or if on the last line
            I3(length(I3)+1)=k-m-(len/2)+2*m; % goes up one line
        end
        if row==1 && mod(m,2)==0
             I3(length(I3)+1)=k+(len/2-m) + m/2; %goes down one line
        end            
    end
end
    
B=Dt';
B(I3)=ones(size(I3));
Dt=B';





  
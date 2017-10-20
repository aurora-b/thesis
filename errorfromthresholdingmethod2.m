%using hyperbolic tan as our localized fuction
x = -4:0.1:4; 
len=length(x)
y=tanh(x)
%plot(x,y);grid on
%%
%Lift the db3 wavelet to transform it into a second generation wavelet
lsdb3 = liftwave('db3'); %lift wavelet
displs(lsdb3); %see what the lifting scheme yields
%%
%find coefficients
lev   = 2;
wname = lsdb3;
yDec = lwt(y,wname,lev)

ca1=lwtcoef('ca',yDec,lsdb3,2,1)
ca2=lwtcoef('ca',yDec,lsdb3,2,2)
cd1=lwtcoef('cd',yDec,lsdb3,2,1)
cd2=lwtcoef('cd',yDec,lsdb3,2,2)

e=1
I = find(abs(cd1)<e);
cd1(I) = zeros(size(I));
I2 = find(abs(cd2)<e);
cd2(I2) = zeros(size(I2));

% yREC=ilwt(ca2,cd2,lsdb3,2)

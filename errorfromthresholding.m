%using hyperbolic tan as our localized fuction
x = -4:0.1:4; 
len=length(x)
y=tanh(x)
%plot(x,y);grid on
%%
%perform a 2nd gen wavelet transform (in physical space)
lsdb3 = liftwave('db3'); %lift wavelet
% Perform LWT at level 'lev'
lev=3
[cA1, cD1]=lwt(y,lsdb3,1)
[cA2, cD2]=lwt(y,lsdb3,2)
[cA3, cD3] =lwt(y,lsdb3,3)

%%
%perform thresholding
e=1E-2
I = find(abs(cD1)<e);
cD1(I) = zeros(size(I));
I2=find(abs(cD2)<e);
cD2(I2)=zeros(size(I2));
I3=find(abs(cD3)<e);
cD3(I3)=zeros(size(I3));

%%
%reconstruct
yRec1=ilwt(cA1,cD1,lsdb3)
yRec2=ilwt(cA2,cD2,lsdb3,2)
yRec3=ilwt(cA3,cD3,lsdb3,3)
plot(yRec1) 
hold on 
plot(yRec2)
plot(yRec3)
grid on

err1=max(max(abs(y-yRec1))) %Why are the errors the same
err2=max(max(abs(y-yRec2)))
err3=max(max(abs(y-yRec3)))






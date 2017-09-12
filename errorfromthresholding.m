%using hyperbolic tan as our localized fuction
x = -4:0.1:4; 
len=length(x)
y=tanh(x)
plot(x,y);grid on
%%
%perform a 2nd gen wavelet transform (in physical space)
lsdb3 = liftwave('db3'); %lift wave
displs(lsdb3); %see what the lifting scheme yields
%use lwtcoef https://www.mathworks.com/help/wavelet/ref/lwtcoef.html
% Perform LWT at level n
lev=3
yDec = lwt(y,lsdb3,lev)
%extract coefficients
for k = 1:lev
        d = lwtcoef('d',yDec,lsdb3,lev,k);
        d = d(:)';
        d = d(ones(1,2^k),:);
        cfd(k,:) = wkeep1(d(:)',len);
end
cfd=cfd(:);
% d1 = abs(lwtcoef('d',yDec,lsdb3,3,1))
% d2= abs(lwtcoef('d',yDec,lsdb3,3,2))
% d3= abs(lwtcoef('d',yDec,lsdb3,3,3))
%%
%perform thresholding
e=1E-3

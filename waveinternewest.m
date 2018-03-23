function [A,D]= waveinternewest(y)
if mod(length(y)-1,2)~=0
    msg = 'Error occurred. Length of initial function -1 must be even.';
    error(msg)
end
%A is the approximation, D is the detail, a3 is the coefficient of lifting
%to preserve mean
%m=order of interpolation
d=zeros(1,length(y));
a=zeros(1,length(y));
yn=y(2:end)
%forward transform
for i=1:length(yn)
if (mod(i,2)==0) && (i~=length(yn));%if even and not an end point
    d(i)=0.5*(yn(i)-0.5*(yn(i-1)+yn(i+1))); %define wavelet definition
end
if mod(i,2)==0 && i==length(yn); %if even and is an end point
    d(i)=0.5*(yn(i)-0.5*(yn(i-1)+y(1)));
end
end 
%define approximation
for i=1:length(y)
if mod(i,2)~=0 && i~=1; %if odd and not an end point
    a(i)=y(i)+0.5*(d(i-1)+d(i+1));
end
if mod(i,2)~=0 && i==1
    a(i)=y(i)+0.5*(d(i+1)+d(length(y)));
end
end
A=a(1:2:end);
D=d(2:2:end);
end













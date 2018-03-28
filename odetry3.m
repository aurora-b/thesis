%Associated with function file 'ode45try2'
clear
tend=0.9;
g=9;
n=2^g; %grid points
b=2*pi; %length of x axis
delx= b/n; %width of space stept

delt=0.1*delx;
%visc=delx^2/8
visc=delx^1.2


x= 0:delx:b-delx; %adds delx each time and specifies grid points
uinit=zeros(1,n); %preallocating u
%%
for i=1:n
    uinit(i)= sin(x(i));
end
%%
[t,u]=ode45(@ode45try2,0:delt:tend,uinit);

w=length(u(:,1))
% 
% final = [u(j,:)]
% middle = [u(j-10,:)]

%%
%Create a movie
for time=1:w
    plot(x,u(time,:));
    title({['1-D Burgers'' equation (\nu = ',num2str(visc),')'];['time(\itt) = ',num2str(delt*time)]});
    grid on;
axis([0,b,-1.5,1.5]);pause(.1)
mov(time)=getframe(figure(1));
end


% v = VideoWriter('burgerrk4.avi');
% v.FrameRate = 100;  % Default 30
% v.Quality = 100;    % Default 75
% open(v)
% writeVideo(v,mov)
% close(v)
%%
%Display just a final plot
plot(x/pi,u(w,:))
grid on;
axis([0 2 -1.5 1.5])
xlabel('x/pi')
ylabel('u(x)')
%%
%plot with x=pi line
h = plot(x,u(w,:)); grid on; axis([0 2*pi -1.5 1.5])
hold on
piy=[-1.5 1.5]
pix= [pi pi]
plot(pix,piy)
xlabel('x')
ylabel('u(x)')



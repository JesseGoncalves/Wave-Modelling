%% Exercise 10.4
% (a)
% compute u on a mesh of points (x,t)
x=-10:0.1:10; t=0:0.01:6; 
[X,T]=meshgrid(x,t);
u=.5*(heaviside(X+2*T).*heaviside(1-X-2*T)+heaviside(X-2*T).*heaviside(1-X+2*T));

% visualize the function through animation
figure; M=moviein(length(t));
for j=1:length(t);
    plot(x,u(j,:)); M(:,j)=getframe;
end;
movie(M)

% (b)
% surface plot
figure; surf(x,t,u)
shading interp 
view(2)

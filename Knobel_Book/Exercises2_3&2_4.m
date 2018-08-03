%% Knobel 2.2: Visualizing 2d Waves
%% Exercise 2.3

% compute u on a mesh of points (x,t)
x=-10:0.1:10; t=0:0.3:6; 
[X,T]=meshgrid(x,t);
u=exp(-(X-T).^2);

% visualize the function through animation
figure; M=moviein(length(t));
for j=1:length(t);
    plot(x,u(j,:)); M(:,j)=getframe;
end;
movie(M)

% slice plot
figure; waterfall(x,t,u)

% surface plot
figure; surf(x,t,u)

% and xt-diagram
figure; pcolor(x,t,u)

%% Exercise 2.4

% compute u on a mesh of points (x,t)
x=-10:0.1:10; t=0:0.1:5;
[X,T]=meshgrid(x,t);
u=exp(-(X-T).^2)+exp(-(X+T).^2);

% visualize with animation
figure; M=moviein(length(t));
for j=1:length(t);
    plot(x,u(j,:)); M(:,j)=getframe;
end;
movie(M)

% slice plot
figure; waterfall(x,t,u)

% surface plot
figure; surf(x,t,u)

% and xt-diagram
figure; pcolor(x,t,u)

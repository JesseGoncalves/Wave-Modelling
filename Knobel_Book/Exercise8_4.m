%% Knobel 8.2: d'Alembert solution of wave equation
% Exercise 8.4
% compute u on a mesh of points (x,t)
x=0:0.3:10; t=0:0.03:10; 
[X,T]=meshgrid(x,t);
u=cos(T).*sin(X);

% visualize the function through animation
figure; M=moviein(length(t));
for j=1:length(t);
    plot(x,u(j,:)); axis([0,10,-1,1]); M(:,j)=getframe;
end;
movie(M)


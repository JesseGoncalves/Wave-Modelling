%% Exercise 4.6
c=0.1; % speed
x=-10:0.1:10; t=0:0.5:20; 
[X,T]=meshgrid(x,t); % spacetime grid
u=4*atan(exp((X-c*T)/sqrt(1-c^2))); % wave equation
M=moviein(length(t)); % create animation
for j=1:length(t);
    plot(x,u(j,:)); M(:,j)=getframe;
end;
movie(M)

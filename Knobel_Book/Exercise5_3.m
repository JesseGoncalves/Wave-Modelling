%% Exercise 5.3
x=-10:0.1:10; t=-10:0.1:10;
[X,T]=meshgrid(x,t);
k1=1; k2=2;
u1=exp(k1^3*T-k1*X);
u2=exp(k2^3*T-k2*X);
A=(k1-k2)^2/(k1+k2)^2;
u=12*(k1^2*u1+k2^2*u2+2*(k1-k2)^2*u1*u2+A*u1*u2*(k1^2*u2+k2^2*u1))/...
   (1+u1+u2+A*u1*u2)^2; % double soliton solution to KdV equation

% animation
figure; M=moviein(length(t));
for j=1:length(t);
    plot(x,u(j,:)); M(:,j)=getframe;
end;
movie(M)

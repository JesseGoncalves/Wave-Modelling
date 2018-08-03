function [] = prob_13_7_d()
% 
% Written by:
% -- 
% John L. Weatherwax                2005-08-14
% 
% email: wax@alum.mit.edu
% 
% Please send comments and especially bug reports to the
% above email address.
% 
%-----

% Holds the states (v_star, u_star): 
q_star = [ 1 1 ]; 

v = linspace( -3, 5, 1000 ); 

% Part (i):
%
up = q_star(2) + sqrt( - ( (-exp(v)) - (-exp(q_star(1))) ) ./ ( v - q_star(1) ) ) .* ( v - q_star(1) ); 
um = q_star(2) - sqrt( - ( (-exp(v)) - (-exp(q_star(1))) ) ./ ( v - q_star(1) ) ) .* ( v - q_star(1) ); 

figure; hp = plot( v, up, '-bo' ); hold on; 
hm = plot( v, um, '-rx' ); grid on; 
plot( q_star(1), q_star(2), 'ok', 'MarkerSize', 10, 'MarkerFaceColor', 'k' ); 
xlabel( 'v' ); ylabel( 'u' ); 
legend( [ hp, hm ], 'Plus Sign', 'Minus Sign' ); 
saveas( gcf, 'chap_13_prob_7_pt_d_i.eps', 'psc2' ); 

% Part (ii):
%
p_star = -(2*q_star(1)+0.1*exp(q_star(1)));
p = -(2*v+0.1*exp(v));
up = q_star(2) + sqrt( - ( p - p_star ) ./ (v - q_star(1) ) ) .* ( v - q_star(1) ); 
um = q_star(2) - sqrt( - ( p - p_star ) ./ (v - q_star(1) ) ) .* ( v - q_star(1) ); 

figure; hp = plot( v, up, '-bo' ); hold on; 
hm = plot( v, um, '-rx' ); grid on; 
plot( q_star(1), q_star(2), 'ok', 'MarkerSize', 10, 'MarkerFaceColor', 'k' ); 
xlabel( 'v' ); ylabel( 'u' ); 
legend( [ hp, hm ], 'Plus Sign', 'Minus Sign' ); 
saveas( gcf, 'chap_13_prob_7_pt_d_ii.eps', 'psc2' ); 

% Part (iii):
%
p_star = -(2*q_star(1));
p = -(2*v);
up = q_star(2) + sqrt( - ( p - p_star ) ./ (v - q_star(1) ) ) .* ( v - q_star(1) ); 
um = q_star(2) - sqrt( - ( p - p_star ) ./ (v - q_star(1) ) ) .* ( v - q_star(1) ); 

figure; hp = plot( v, up, '-bo' ); hold on; 
hm = plot( v, um, '-rx' ); grid on; 
plot( q_star(1), q_star(2), 'ok', 'MarkerSize', 10, 'MarkerFaceColor', 'k' ); 
xlabel( 'v' ); ylabel( 'u' ); 
legend( [ hp, hm ], 'Plus Sign', 'Minus Sign' ); 
saveas( gcf, 'chap_13_prob_7_pt_d_iii.eps', 'psc2' ); 

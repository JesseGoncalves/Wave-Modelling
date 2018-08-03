function [] = prob_13_7_e()
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

% Compute the {\em left} Rankine-Hugoniot curve (valid for v_m > v_l):
ql_star = [ 1 1 ]; 
v = linspace( -3, 5, 200 ); 
um = ql_star(2) - sqrt( - ( (-exp(v)) - (-exp(ql_star(1))) ) ./ ( v - ql_star(1) ) ) .* ( v - ql_star(1) ); 

% Compute the {\em right} Rankine-Hugoniot curve (valid for v_m > v_r):
qr_star = [ 3 4 ]; 
v = linspace( -3, 5, 200 ); 
up = qr_star(2) + sqrt( - ( (-exp(v)) - (-exp(qr_star(1))) ) ./ ( v - qr_star(1) ) ) .* ( v - qr_star(1) ); 

% Plot them: 
figure; 

hm = plot( v, um, '-rx' ); hold on; 
plot( ql_star(1), ql_star(2), 'ok', 'MarkerSize', 10, 'MarkerFaceColor', 'k' ); 

hp = plot( v, up, '-bo' ); grid on; 
plot( qr_star(1), qr_star(2), 'ok', 'MarkerSize', 10, 'MarkerFaceColor', 'k' ); 

xlabel( 'v' ); ylabel( 'u' ); title( 'Two-Shock Solution (Maynot satistisfy the entropy condition)' ); 
legend( [ hp, hm ], '2-shock curver', '1-shock curve' );

fprintf('Looks like they intersect at (v,u) = (%10.6f, %10.6f)\n',1.7,-0.4); 
hold on; plot( 1.7, -0.4, 'dk', 'MarkerSize', 8, 'MarkerFaceColor', 'k' ); 
%saveas( gcf, 'chap_13_prob_7_pt_e.eps', 'psc2' ); 

% Solving with Matlabs fzero: 
%
v_m = fzero( 'prob_13_7_e_fn', 1.7 );
u_m = ql_star(2) - sqrt( - ( (-exp(v_m)) - (-exp(ql_star(1))) ) ./ ( v_m - ql_star(1) ) ) .* ( v_m - ql_star(1) ); 
fprintf('Actually they intersect at (v,u) = (%10.6f, %10.6f)\n',v_m,u_m);
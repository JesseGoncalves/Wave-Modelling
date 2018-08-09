function [u_l_new,u_r_new] = burgersRiemannFVM(u_l,u_r,dx,dt);
% solves u_t + uu_x = 0 with initial data u_l for x<0 and u_r for x>0
% dx = cell-size for FVM
% dt = time-step for FVM


if u_l == u_r
    u_l_new = u_l; u_r_new = u_r;
elseif u_l > u_r % shock wave
    s = (u_r^2 - u_l^2)/(2*(u_r - u_l)); % shock speed
    if dx < s*dt || dx > abs(s*dt)
        u_l_new = u_l;
    else
        u_l_new = u_r;
    end
    if dx < s*dt
        u_r_new = u_l;
    else
        u_r_new = u_r;
    end
else % rarefaction wave
    if dx < u_l*dt || abs(u_l*dt) < dx
        u_l_new = u_l;
    else
        u_l_new = -dx/dt;
    end
    if dx > u_r*dt
        u_r_new = u_r;
    else
        u_r_new = dx/dt;
    end
end
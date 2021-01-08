% Writer: Nugraha Setya Ardi
function [delta, rp, e_m, e_p, a_m, a_p, delta_vp, vp_m, vp_p] = flybyPow(v_infM, v_infP, muP)
    delta = acos(dot(v_infM,v_infP)/(norm(v_infM)*norm(v_infP)));
    
    % Setting initial guess for solving rp
    e = 1/sin(delta/2);
    a1 = -muP/norm(v_infP)^2;
    a2 = -muP/norm(v_infM)^2;   
    % Solving rp
    fun = @(rp) asin(1/(1+rp*norm(v_infM)^2/muP))+asin(1/(1+rp*norm(v_infP)^2/muP))-delta;
    rp_guess1 = a1*(1-e);
    rp_guess2 = a2*(1-e);
    rp_guess = (rp_guess1 + rp_guess2)/2;
    rp = fzero(fun,rp_guess);
    
    e_m = 1 + rp*norm(v_infM)^2/muP;
    e_p = 1 + rp*norm(v_infP)^2/muP;
    a_m = rp/(1-e_m);
    a_p = rp/(1-e_p);
    vp_m = sqrt(norm(v_infM)^2 + 2*muP/rp);
    vp_p = sqrt(norm(v_infP)^2 + 2*muP/rp);
    delta_vp = abs(norm(vp_p) - norm(vp_m));
end

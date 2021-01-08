function [delta, rp, e_m, e_p, a_m, a_p, delta_vp, vp_m, vp_p] = flybyPow(v_infM, v_infP, muP)
    % Writer: Nugraha Setya Ardi
    % This function calculates various parameters during powered gravity assist/fly-by (which means that some amounts of delta velocity is given to 
    % the spacecraft when its position is in hyperbolic pericenter)
    % Input: 
    %   1. v_infM : spacecraft's infinite velocity vector in planetocentric before fly-by. [km/s]
    %   2. v_infP : spacecraft's infinite velocity vector in planetocentric after fly-by. [km/s]
    %   3. muP : standar gravitational parameter of planet fly-by. [km^3/s^2]
    % Output:
    %   1. delta : deflection angles. [rad]
    %   2. rp : pericenter of hyperbolic trajectory [km]
    %   3. e_m : eccentricity of hyperbolic before reaching pericenter.
    %   4. e_p : eccentricity of hyperbolic after reaching pericenter.
    %   5. a_m : semi-major axis of hyperbolic before reaching pericenter. [km]
    %   6. a_p : semi-major axis of hyperbolic before reaching pericenter. [km]
    %   7. delta_vp : delta velocity required to do powered fly-by. [km/s]
    %   8. vp_m : velocity at pericenter before being given delta_vp. [km/s]
    %   9. vp_p : velocity at pericenter after being given delta_vp. [km/s]
    
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

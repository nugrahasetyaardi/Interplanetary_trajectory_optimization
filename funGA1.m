%% Fitness Function for Genetic Algorithm and Particle Swarm Optimization
% Writer: Nugraha Setya Ardi
function dVtot = funGA1(time_window)
    muSun = astroConstants(4);
    muP = astroConstants(14);
    t_dep = time_window(1);
    t_ga = time_window(2);
    t_arr = time_window(3);
    Dt1 = (t_ga-t_dep)*24*60*60;
    Dt2 = (t_arr-t_ga)*24*60*60;
    
    if Dt1 <= 1*365*24*3600 || Dt2 <= 1*365*24*3600 % Condition for negative time of flight
        dVtot = 1996;
    else
        [kep1,kSun] = uplanet(t_dep,1);
        [kep2,kSun] = uplanet(t_ga,4);
        [kep3,kSun] = uplanet(t_arr,5);
        [r1,v1] = kep2car(kep1,muSun);
        [r2,v2] = kep2car(kep2,muSun);
        [r3,v3] = kep2car(kep3,muSun);
        [A,P,E,ERROR,v1l,v2l,TPAR,THETA] = lambertMR( r1, r2, Dt1, muSun, 0, 0, 1 );
        [A,P,E,ERROR,v2ll,v3l,TPAR,THETA] = lambertMR( r2, r3, Dt2, muSun, 0, 0, 1 );
        v_infM = v2l - v2;
        v_infP = v2ll - v2;
        [delta, rp, e_m, e_p, a_m, a_p, delta_vp, vp_m, vp_p] = flybyPow(v_infM, v_infP, muP);
        deltaV_dep = norm(v1-v1l);
        deltaV_ga = delta_vp;
        deltaV_arr = norm(v3-v3l);
        dVtot = deltaV_dep + deltaV_arr + deltaV_ga;
            if rp < astroConstants(24)+ 300; %mars radius + its atmosphere
                dVtot = 1000;
            end
    end
end



    
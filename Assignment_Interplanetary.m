%% Interplanetary Gravity Assist
clear all

% Departure: Mercury
% Flyby: Mars
% Arrival: Jupiter
% Earliest departure: 2025/06/01
% Latest arrival: 2065/06/01

muSun = astroConstants(4);
muP = astroConstants(14);

tdep_window = [2025,6,1,0,0,0;2030,6,1,0,0,0];
tga_window = [2035,6,1,0,0,0;2045,6,1,0,0,0];
tarr_window = [2050,6,1,0,0,0;2065,6,1,0,0,0];


t1_start_mjd2000 = date2mjd2000(tdep_window(1,1:6));
t1_end_mjd2000 = date2mjd2000(tdep_window(2,1:6));
t2_start_mjd2000 = date2mjd2000(tga_window(1,1:6));
t2_end_mjd2000 = date2mjd2000(tga_window(2,1:6));
t3_start_mjd2000 = date2mjd2000(tarr_window(1,1:6));
t3_end_mjd2000 = date2mjd2000(tarr_window(2,1:6));

% Iteration
ndep = t1_end_mjd2000 - t1_start_mjd2000;
nga = t2_end_mjd2000 - t2_start_mjd2000;
narr = t3_end_mjd2000 - t3_start_mjd2000;
t_dep = linspace(t1_start_mjd2000,t1_end_mjd2000,10);
t_ga = linspace(t2_start_mjd2000,t2_end_mjd2000,10);
t_arr = linspace(t3_start_mjd2000,t3_end_mjd2000,10);

min_dV = 1000;
for i =1:10
    for j=1:10
        for k=1:10
            Dt1 = (t_ga(j)-t_dep(i))*24*60*60;
            Dt2 = (t_arr(k)-t_ga(j))*24*60*60;
            a1(i,j,k) = Dt1;
            a2(i,j,k) = Dt2;
            [kep1,kSun] = uplanet(t_dep(i),1);
            [kep2,kSun] = uplanet(t_ga(j),4);
            [kep3,kSun] = uplanet(t_arr(k),5);
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
            deltaV_tot(i,j,k) = deltaV_dep + deltaV_arr + deltaV_ga;
            
            if min_dV >= deltaV_dep + deltaV_arr + deltaV_ga
                min_dV = deltaV_dep + deltaV_arr + deltaV_ga;
                date_dep = mjd20002date(t_dep(i));
                date_ga = mjd20002date(t_ga(j));
                date_arr = mjd20002date(t_arr(k));
            end
            
        end
    end
end



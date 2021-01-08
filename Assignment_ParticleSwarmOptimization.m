% Interplanetary Gravity Assist
% Contributor: Nugraha Setya Ardi
clear all

% Departure: Mercury   
% Flyby: Mars
% Arrival: Jupiter
% Earliest departure: 2025/06/01
% Latest arrival: 2065/06/01

muSun = astroConstants(4);
muP = astroConstants(14);

tdep_window = [2025,6,1,0,0,0;2065,6,1,0,0,0];
tga_window = [2025,6,1,0,0,0;2065,6,1,0,0,0];
tarr_window = [2025,6,1,0,0,0;2065,6,1,0,0,0];

t1_start_mjd2000 = date2mjd2000(tdep_window(1,1:6));
t1_end_mjd2000 = date2mjd2000(tdep_window(2,1:6));
t2_start_mjd2000 = date2mjd2000(tga_window(1,1:6));
t2_end_mjd2000 = date2mjd2000(tga_window(2,1:6));
t3_start_mjd2000 = date2mjd2000(tarr_window(1,1:6));
t3_end_mjd2000 = date2mjd2000(tarr_window(2,1:6));

fitness = @funGA1;
nvars = 3;
lb = [t1_start_mjd2000 t2_start_mjd2000 t3_start_mjd2000]
ub = [t1_end_mjd2000 t2_end_mjd2000 t3_end_mjd2000]

%options = optimoptions('particleswarm','SwarmSize',5000);
options = optimoptions('particleswarm','SwarmSize',10000)%'HybridFcn',@fmincon);
%[x,fval,exitflag] = particleswarm(fun,nvars,lb,ub,options)

[time_window, dVtot] = particleswarm(fitness,nvars,lb,ub,options);
date_dep = mjd20002date(time_window(1))
date_flyby = mjd20002date(time_window(2))
date_arr = mjd20002date(time_window(3))

%% Orbit propagation on Heliocentric
Dt1 = time_window(2) - time_window(1);
Dt2 = time_window(3) - time_window(2);
[kep1,kSun] = uplanet(time_window(1),1);
[kep2,kSun] = uplanet(time_window(2),4);
[kep3,kSun] = uplanet(time_window(3),5);
[r1,v1] = kep2car(kep1,muSun);
[r2,v2] = kep2car(kep2,muSun);
[r3,v3] = kep2car(kep3,muSun);

[A,P,E,ERROR,v1l,v2l,TPAR,THETA] = lambertMR( r1, r2, Dt1*24*3600, muSun, 0, 0, 1 );
[A,P,E,ERROR,v2ll,v3l,TPAR,THETA] = lambertMR( r2, r3, Dt2*24*3600, muSun, 0, 0, 1 )

tRange1 = linspace(0,Dt1*24*3600,10000);
F1 = [r1.';v1l.'];
[tSol1,FSol1] = ode45( @(t,F1) SPFun(F1,muSun), tRange1,F1);

tRange2 = linspace(0,Dt2*24*3600,1000);
F2 = [r2.';v2ll.'];
[tSol2,FSol2] = ode45( @(t,F2) SPFun(F2,muSun), tRange2,F2);

figure
plot3(r2(1),r2(2),r2(3),'O')
hold on
plot3(FSol1(:,1),FSol1(:,2),FSol1(:,3))
hold on
plot3(FSol2(:,1),FSol2(:,2),FSol2(:,3),'black')
plot3(0,0,0,'O')
hold on
plot3(r2(1),r2(2),r2(3),'O')
hold on
plot3(r1(1),r1(2),r1(3),'O')
hold on
plot3(r3(1),r3(2),r3(3),'O')
text(0,0,0,'Sun')
text(r1(1),r1(2),r1(3),'Merc @ dep')
text(r2(1),r2(2),r2(3),'Mars @ flyby')
text(r3(1),r3(2),r3(3),'Jup @ arr')
%% Checking....

t_dep1 = date2mjd2000([2035,4,8,5,1,5.019334853]);
t_ga1 = date2mjd2000([2039,2,10,14,14,27.2829777]);
t_arr1 = date2mjd2000([2041,10,20,8,42,50.42723]);
time_window1 = [t_dep1,t_ga1,t_arr1];
dVtot = funGA1(time_window1)


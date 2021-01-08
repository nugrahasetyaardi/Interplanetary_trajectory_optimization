function [r,v] = kep2car(coe,mu)
a = coe(1);
e = coe(2);
i = coe(3);
RA = coe(4);
w = coe(5);
f = coe(6);

h=sqrt(mu*a*(1-e^2));
rp = (h^2/mu) * (1/(1 + e*cos(f))) * (cos(f)*[1;0;0] + sin(f)*[0;1;0]);
vp = (mu/h) * (-sin(f)*[1;0;0] + (e + cos(f))*[0;1;0]);

R3_W = [ cos(RA) sin(RA) 0; -sin(RA) cos(RA) 0; 0 0 1];
R1_i= [1 0 0; 0 cos(i) sin(i); 0 -sin(i) cos(i)];
R3_w = [ cos(w) sin(w) 0; -sin(w) cos(w) 0; 0 0 1];

Q_pX = (R3_w*R1_i*R3_W);

r = Q_pX'*rp;
v = Q_pX'*vp;
% Convert r and v into row vectors:
r = r';
v = v';
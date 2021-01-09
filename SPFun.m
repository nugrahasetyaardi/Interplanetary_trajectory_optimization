function dFdt = SPFun(F,mu)
    % Writer: Nugraha Setya Ardi
    % This function is a function for integrating numerically newton's gravity equation between 2 celestial bodies, assuming m1 >> m2 (so baricenter lies at the center of m1).
    % You can specify disturbance acceleration in adx, ady, and adz
    x = F(1);
    y = F(2);
    z = F(3);
    vx = F(4);
    vy = F(5);
    vz = F(6);
    
    rR=sqrt(x^2+y^2+z^2);
    adx = 0;
    ady = 0;
    adz = 0;
    
    dxdt = vx;
    dydt = vy;
    dzdt = vz;
    dvxdt = (-mu/(rR^3))*x+adx;
    dvydt = (-mu/(rR^3))*y+ady;
    dvzdt = (-mu/(rR^3))*z+adz;
    dFdt = [dxdt;dydt;dzdt;dvxdt;dvydt;dvzdt];

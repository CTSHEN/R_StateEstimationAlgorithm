%%
% Convert ecef position to Latitude, Longitude and Altitude
% Using Ferrari's method

function LLA = POS_ECEF2LLA(pos)
    X = pos(1);
    Y = pos(2);
    Z = pos(3);

    a=6378137; % Earth equatorial radius
    a_sq = a^2;
    e = 8.181919084261345e-2;
    e_sq = 6.69437999014e-3;

    f = 1/298.257223563;
    b = a*(1-f);

    r = sqrt(X^2 + Y^2);
    ep_sq = (a^2 - b^2)/b^2;
    ee = (a^2 - b^2);
    F = (54*b^2)*Z^2;
    g = r^2 + (1-e_sq)*Z^2 - e_sq*ee*2;
    c = e_sq^2*F*r^2/g^3;
    s = (1+c+sqrt(c^2 + 2*c))^(1/3);
    k = s + 1 + 1/s;
    p = F/(3*k^2*g^2);
    q = sqrt(1 + 2*p*e_sq^2);
    r_0 = -(p*e_sq*r)/(1+q) + sqrt(0.5*a^2*(1+(1/q)) - p*Z^2*(1 - e_sq)/(q*(q+q)) - 0.5*p*r^2);
    u = sqrt((r - e_sq*r_0)^2 + Z^2);
    v = sqrt((r - e_sq*r_0)^2 + (1 - e_sq)*Z^2);
    z_0 = b^2*Z/(a*v);
    h = u*(1 - b^2/(a*v));

    phi = atan((Z + ep_sq*z_0)/r);
    lambd = atan2(Y, X);

    % return LLA in degree, degree, meters
    LLA = [ rad2deg(phi), rad2deg(lambd), h];
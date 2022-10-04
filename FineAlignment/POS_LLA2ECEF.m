%%
% Return a 3x1 vector ofposition vector in ECEF frame with Latitude,
% Longitude and Altitude inputs.
function P_ECEF = POS_LLA2ECEF(lat, long ,alt)
f = 1/298.257223563;
R = 6378137;

ls = atan((1-f)^2*tand(lat));
rs =sqrt(R^2/(1+(1/(1-f)^2 - 1)*sin(ls)^2));

P_ECEF = [ rs*cos(ls)*cosd(long)+alt*cosd(lat)*cosd(long);
                  rs*cos(ls)*sind(long)+alt*cosd(lat)*sind(long);
                  rs*sin(ls)+alt*sind(lat)];
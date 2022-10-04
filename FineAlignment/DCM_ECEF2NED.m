function DCM =DCM_ECEF2NED(Lat, Long)
cLat = cosd(Lat);
cLong = cosd(Long);
sLat = sind(Lat);
sLong = sind(Long);

DCM = [-sLat*cLong  -sLat*sLong  cLat;
              -sLong  cLong  0;
              -cLat*cLong  -cLat*sLong  -sLat];
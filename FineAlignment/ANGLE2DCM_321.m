%%
% This function calculates the dcm matrix from local NED frame to body
% frame with yaw-pitch-roll rotation sequence.

function DCM = ANGLE2DCM_321(yaw, pitch, roll)
c1 = cos(yaw);
c2 = cos(pitch);
c3 = cos(roll);
s1 = sin(yaw);
s2 = sin(pitch);
s3 = sin(roll);

DCM = [c1*c2  c1*s2*s3-c3*s1  c1*s2*c3+s3*s1;
              s1*c2  s1*s2*s3+c3*c1  s1*s2*c3-s3*c1;
              -s2  c2*s3  c2*c3]';


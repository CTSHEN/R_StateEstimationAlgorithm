% Static TRIAD

% Define some constants
EarthRotation = 7.2921150e-5; % Earth rotation angular speed, in rad/sec.
e = 0.0818191908425;  % eccentricity of the ellipsoid.

TestLat = 24.72534; %24 + 43/60 + 29/3600;  % Latitude of TiSpace  29
TestLong = 120.91023; %120 + 54/60 + 35/3600;  % Longitude of TiSpace  35
g = 9.7803253359*(1+0.001931853*sin(deg2rad(TestLat))^2)/sqrt(1-e^2*sin(deg2rad(TestLat))^2);  %m/s^2

KVH_ARW = 0.012;  % Angle Rnadon Walk, in deg/sqrt(hr).
KVH_VRW = 0.007;  % Velocity Random Walk, in mg/sqrt(Hz).

% Deal with units of ARW and VRW
KVH_ARW_newUnit = KVH_ARW*pi()/(180 * 3600);  % the unit is rad/sqrt(sec).
KVH_VRW_newUnit = KVH_VRW * 1e-6 * g;  % the unit is (m/s^2)*(1/sqrt(Hz)).

% The refernece vectors
w_IE_E = [0;0;EarthRotation];  % Since ECI and ECEF frames are coincidance in 3rd axis and
                                                     %  w_IE has only one
                                                     %  value on 3rd
                                                     %  element, therefore
                                                     %  w_IE_E = w_IE_I.

E2LD = dcmecef2ned(TestLat, TestLong);  % A DCM from ecef to local NED frame
LD2LU = [0 1 0
                 1 0 0
                 0 0 -1];   % DCM from NED to ENU.


g_LD = [0;0;-g];  % Reference gravitational measurement in local NED frame.
% CHECK THE DIRECTION OF GRAVITY

g_E = E2LD' * g_LD;  % Transform g from NED to ECEF.
g_LU = LD2LU*g_LD;
w_IE_LU = LD2LU*E2LD*w_IE_E;
w_IE_LD = E2LD*w_IE_E;

r1 = w_IE_LU;  % reference vector 1

r2 = g_LU;  % reference vector 2
%% Long time stationary case
Start_idx = 33764;%200*60*10 + 1; 
End_idx =55914;%200*60*60;
if End_idx > size(accM_enu,2)
    End_idx = size(accM_enu,2);
end

DataSamples = End_idx - Start_idx;

b1 = [mean(rateM_enu(1,Start_idx:End_idx)); mean(rateM_enu(2,Start_idx:End_idx)); mean(rateM_enu(3,Start_idx:End_idx)) ];
b2 = [mean(accM_enu(1,Start_idx:End_idx)); mean(accM_enu(2,Start_idx:End_idx)); mean(accM_enu(3,Start_idx:End_idx)) ];

%% remove bias
% b1 =[-3.93E-05;-5.34E-05;3.07E-05];
% b2 = [0; 0; 9.81];


% gbias = [3.85E-06;-9.22E-07;-9.76E-07];  % average bias of first 7 exps 30C
% abias = [0.01877;-0.01173;0.00778];  % average bias of first 7 exps 30C

% gbias = [-3.54E-06;1.57E-06;8.00E-06];  % average bias of first 5 exps 39C
% abias = [0.019023033;-0.014227794;0.008211323 ]; % average bias of first 5 exps 39C

gbias = [-9.507704E-07;-3.886855E-06;2.60115E-06];  % interpolate bias test
abias = [0.016612834;-0.011980425;0.006255176]; % interpolate bias test

% abias = [0; 0; 0];
% gbias = [0; 0; 0];
b1_b = b1 - gbias;
b2_b = b2 - abias;


%% TRIAD
r1 = r1/norm(r1);
r2 = r2/norm(r2);
 b1_b = b1_b/norm(b1_b);
b2_b = b2_b/norm(b2_b);

bx = cross(b1_b, b2_b);
bx = bx / norm(bx);
rx =  cross(r1, r2);
rx = rx/norm(rx);
bxx = cross(b2_b,bx);
bxx = bxx/norm(bxx);
rxx = cross(r2,rx);
rxx = rxx/norm(rxx);

    BB = [bxx b2_b bx];
    RR = [rxx r2 rx];
%     BB = [b1 b2 bx];
%     RR = [r1 r2 rx];
    AA_LU = BB*RR';
    [AALU1, AALU2, AALU3] = dcm2angle(AA_LU,'XYZ');

    AALU3_deg = rad2deg(AALU3);
    AALU2_deg = rad2deg(AALU2);
    AALU1_deg = rad2deg(AALU1);

    AA_LD = AA_LU*LD2LU;
    [AALD1, AALD2, AALD3] = dcm2angle(AA_LD,'ZYX');

    AALD3_deg = rad2deg(AALD3);
    AALD2_deg = rad2deg(AALD2);
    AALD1_deg = rad2deg(AALD1);

%% Stationary, different azimuth case, TRIAD
for idx = 1:size(accM_enu,2)
    b1 = rateM_enu(:,idx);
    b2 = accM_enu(:,idx);

    r1 = r1/norm(r1);
    r2 = r2/norm(r2);
    b1 = b1/norm(b1);
    b2 = b2/norm(b2);

    bx = cross(b1, b2);
%     bx = bx / norm(bx);
    rx =  cross(r1, r2);
%     rx = rx/norm(rx);
    bxx = cross(b2,bx);
%     bxx = bxx/norm(bxx);
    rxx = cross(r2,rx);
%     rxx = rxx/norm(rxx);

    BB = [bxx b2 bx];
    RR = [rxx r2 rx];
    
    AA_LU = BB*RR';
    [AALU1(idx), AALU2(idx), AALU3(idx)] = dcm2angle(AA_LU,'XYZ');

    AALU3_deg(idx) = rad2deg(AALU3(idx));
    AALU2_deg(idx) = rad2deg(AALU2(idx));
    AALU1_deg(idx) = rad2deg(AALU1(idx));

    % calculate yaw estimate error
    yawError(idx) = AALU3_deg(idx) - (-180+idx);
end
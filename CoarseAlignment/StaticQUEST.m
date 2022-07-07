%  This program uses QUEST algorithm with two observations, earth rotating
%  angular velocity and gravity, to estimate current attitude.
%  QUEST : Quaternion Estimator
%  This QUEST algorithm is speciallized for two observation vectors.
%
%  The algorithm is based on the book:
%  F. Landis Markley, John L. Crassidis: Fundamentals of Spacecraft Attitude
%  Determination and Control
%
%  20220706
%  v0.1 First Implimentation of QUEST
%
%  Some nomenclatures:
%       I : ECI frame
%       E : ECEF frame
%       LU : Local ENU frame
%       LD : Local NED frame
%       B : Body frame ( IMU frame )
%
%       for vectors:
%           w_IE_I : vector "w" of E frame w.r.t. I frame, express on I
%           frame.
%          

% Define some constants
EarthRotation = 7.2921150e-5; % Earth rotation angular speed, in rad/sec.
g = 9.80665; % gravitation acceleration, in m/s^2

TestLat = 24 + 43/60 + 29/3600;  % Latitude of TiSpace
TestLong = 120 + 54/60 + 35/3600;  % Longitude of TiSpace

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
%w_IE_LU = LD2LU*E2LD*w_IE_E;  % Express w_IE in ENU frame.
w_IE_LU = E2LD*w_IE_E;  % Express w_IE in NED frame.

g_LU = [0;0;g];  % gravitational acceleration in local ENU frame.

r1 = w_IE_LU;  % reference vector 1
%r1 = r1/norm(r1);
r2 = g_LU;  % reference vector 2
%r2 = r2/norm(r2);

%% Generate Test Data

yawTestd = 0;  %deg
pitchTestd = 60;
rollTestd = 0;
% generate ZYX sequence DCM
dcmTest = angle2dcm(deg2rad(yawTestd), deg2rad(pitchTestd), deg2rad(rollTestd),'ZYX');
% generate test vector
b1 = dcmTest*r1;
b2 = dcmTest*r2;

% Choose the data point for testing
Dataidx = 1;

%% 
%  Extract a measurement data to test the QUEST algorithm
% assign measurement vectors

%for Dataidx = 1:size(RateX)

   %b1 = [RateX(Dataidx); RateY(Dataidx); RateZ(Dataidx)];
   %b1 = b1/norm(b1);
   %b2 = [AccX(Dataidx); AccY(Dataidx); AccZ(Dataidx)];
   %b2 = b2/norm(b2);

    % Compute the weight vector
    %a1 = 1/KVH_ARW_newUnit^2;
    %a2 = 1/KVH_VRW_newUnit^2;
    a1 = 1/2;
    a2 = 1/2;
 
    % Compute the quaternion
    bx = cross(b1, b2);
    %bx = bx / norm(bx);
    rx =  cross(r1, r2);
    %rx = rx/norm(rx);

    alpha = (1 + dot(bx,rx)) * (a1*dot(b1, r1) + a2*dot(b2, r2)) + dot(cross(bx, rx), (a1*cross(b1,r1) + a2*cross(b2,r2)) );

    beta = dot( (bx + rx), (a1*cross(b1,r1) + a2*cross(b2,r2)));

    gamma = sqrt(alpha^2 + beta^2);

    if alpha >= 0
        prefix = 1/ (2*sqrt(gamma * (gamma + alpha) * (1 + dot(bx, rx))));
        q_vec = [ (gamma + alpha)*cross(bx,rx) + beta*(bx + rx);
                        (gamma + alpha)*(1 + dot(bx,rx)) ];

        q_hat(:,Dataidx) = prefix * q_vec;   % This q_hat is JPL style.
    else
        prefix = 1/(2*sqrt(gamma*(gamma - alpha)*(1 + dot(bx,rx))));
        q_vec = [ beta*cross(bx,rx) + (gamma - alpha)*(bx + rx);
                        beta*(1 + dot(bx,rx)) ];

        q_hat(:,Dataidx) = prefix * q_vec;  % This q_hat is JPL style.
    end


    % Transform the q_hat from JPL style to Hamilton style, which is the style
    % used in MATLAB

    quat_est = quaternion(q_hat(4,Dataidx), q_hat(1,Dataidx), q_hat(2,Dataidx), q_hat(3,Dataidx));
    % Turn the quaternion into ypr angle.
    YPR(Dataidx,:) = eulerd(quat_est, "ZYX", "point");

%end




%  This program uses QUEST algorithm with two observations, earth rotating
%  angular velocity and gravity, to estimate current attitude.
%  QUEST : Quaternion Estimator
%  This QUEST algorithm is speciallized for two observation vectors.
%  rate and acceleration vector input is in ENU.
%
%  The algorithm is based on the book:
%  F. Landis Markley, John L. Crassidis: Fundamentals of Spacecraft Attitude
%  Determination and Control
%
%  20220706
%  v0.1 First Implimentation of QUEST
% 20220707
% v0.11 Set ECEF as the reference frame instead of local ENU frame.
% v0.12 Check some quaternion normalization to avoid numerical problem.
% v0.13 Add estimatino of covariance of error angle.
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

DataSamples = 2000;  % How many data samples used in the test

TestLat = 24.72534; %24 + 43/60 + 29/3600;  % Latitude of TiSpace  29
TestLong = 120.91023; %120 + 54/60 + 35/3600;  % Longitude of TiSpace  35

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

r1 = w_IE_LU;  % reference vector 1, LU
r2 = g_LU;  % reference vector 2 , LU

% r1 = w_IE_E;  % reference vector 1, E
% r2 = g_E;  % reference vector 2, E




%% Generate Test Data

yawTestd_LU = 0;  %deg
pitchTestd_LU =10;
rollTestd_LU = 0;
% generate ZYX sequence DCM
dcm_LU2B = angle2dcm(deg2rad(rollTestd_LU), deg2rad(pitchTestd_LU), deg2rad(yawTestd_LU),'XYZ');
% generate test vector
b1 = dcm_LU2B*r1;  % r1 = w)IE)LU
b2 = dcm_LU2B*LD2LU*g_LD;

% Choose the data point for testing
Dataidx = 1;

%% Using mean value as measurement
% b1 = [mean(RateX(:)); mean(RateY(:)); mean(RateZ(:)) ];
% b2 = [mean(AccX(:)); mean(AccY(:)); mean(AccZ(:)) ];
Start_idx = 200*60*0 + 1; 
End_idx = 200*60*60;
if End_idx > size(accM_enu,2)
    End_idx = size(accM_enu,2);
end

DataSamples = End_idx - Start_idx;

b1 = [mean(rateM_enu(1,Start_idx:End_idx)); mean(rateM_enu(2,Start_idx:End_idx)); mean(rateM_enu(3,Start_idx:End_idx)) ];
b2 = [mean(accM_enu(1,Start_idx:End_idx)); mean(accM_enu(2,Start_idx:End_idx)); mean(accM_enu(3,Start_idx:End_idx)) ];

% GAbias = States(5:10, 720000);
% b1_bias = GAbias(1:3);
% b2_bias = GAbias(4:6);
% b1_bias =  [deg2rad(-2e-4) ; deg2rad(0.9e-4) ; deg2rad(-2.25e-4)];
% b2_bias = [ -0.006; -0.012; 0.01];
% b1 = b1-b1_bias;
% b2 = b2-b2_bias;



%% 
% assign measurement vectors
r1 = r1/norm(r1);
r2 = r2/norm(r2);

for Dataidx = 1:1 %DataSamples

   
   %b1 = [RateX(Dataidx); RateY(Dataidx); RateZ(Dataidx)];
   
   
   b1 = b1/norm(b1);
   %b2 = [AccX(Dataidx); AccY(Dataidx); AccZ(Dataidx)];
   
   b2 = b2/norm(b2);

   % Compute the weight vector
    a1 = 1/KVH_ARW_newUnit^2;
    a2 = 1/KVH_VRW_newUnit^2;
   
    % Compute the quaternion
    bx = cross(b1, b2);
    bx = bx / norm(bx);
    rx =  cross(r1, r2);
    rx = rx/norm(rx);
    bxx = cross(b2,bx);
%     bxx = bxx/norm(bxx);
    rxx = cross(r2,rx);
%     rxx = rxx/norm(rxx);

    alpha = (1 + dot(bx,rx)) * (a1*dot(bxx, rxx) + a2*dot(b2, r2)) + ...
        dot(cross(bx, rx), (a1*cross(bxx,rxx) + a2*cross(b2,r2)) );

    beta = dot( (bx + rx), (a1*cross(bxx,rxx) + a2*cross(b2,r2)));

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

    quat_E2B = [ q_hat(4,Dataidx)  q_hat(1,Dataidx)  q_hat(2,Dataidx)  q_hat(3,Dataidx) ];
    quat_E2B = quat_E2B/norm(quat_E2B);
    
    % In this comment section, "*" denotes quaternion multiplication.
    % quat_E2B = quat_E2LD * quat_LD2LU * quat_LU2B 
    % therefore, quat_LU2B = (quat_E2LD * quat_LD2LU)^(-1) * quat_E2B
    quat_LD2LU = dcm2quat(LD2LU);
    quat_E2LD = dcm2quat(E2LD);

%     quat_LU2B = quatmultiply( quatinv(quatmultiply(quat_E2LD, quat_LD2LU)), quat_E2B);
%     quat_LU2B = quatmultiply( quat_E2B, quatinv(quatmultiply( quat_LD2LU,quat_E2LD)));
%     quat_LU2B = quatmultiply( quatmultiply(  quat_E2B,quatinv(quat_E2LD)),quatinv(quat_LD2LU) );
   
    
    % Turn the quaternion into ypr angle.
    
%     RPY(Dataidx,:) = rad2deg( quat2eul(quat_LU2B, 'XYZ' ));
    RPY(Dataidx,:) = rad2deg( quat2eul([q_hat(4,Dataidx)  q_hat(1,Dataidx)  q_hat(2,Dataidx)  q_hat(3,Dataidx)], 'XYZ' ));

    % Compute error angle covariance
    P_errorAngle = (KVH_VRW_newUnit^2 *b1*b1' + KVH_ARW_newUnit^2 *b2*b2')/ norm(cross(b1,b2))^2 +...
        ((KVH_ARW_newUnit^2 * KVH_VRW_newUnit^2)/(KVH_ARW_newUnit^2 + KVH_VRW_newUnit^2)) * bx * bx';

    % Error angle variance (EAV)
    EAV(Dataidx,:) = sqrt( diag(P_errorAngle)' );

end




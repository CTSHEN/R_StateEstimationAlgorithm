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
w_IE_LU = LD2LU*E2LD*w_IE_E;  % Express w_IE in ENU frame.

g_LU = [0;0;g*1e-5];  % gravitational acceleration in local ENU frame.

r1 = w_IE_LU;  % reference vector 1

r2 = g_LU;  % reference vector 2

%% Generate Test Data

yawTestd = 135;  %deg
pitchTestd = 0;
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
   
   %b2 = [AccX(Dataidx); AccY(Dataidx); AccZ(Dataidx)];
   
   % Compute the weight vector
    sigma1 = 1/KVH_ARW_newUnit^2;
    sigma2 = 1/KVH_VRW_newUnit^2;
    sigma_tot_sq = ((sigma1^2)^(-1) + (sigma2^2)^(-1))^(-1);

    a1 = sigma_tot_sq/sigma1^2;
    a2 = sigma_tot_sq/sigma2^2;

    % compute attitude profile (B) matrix
    B = a1*b1*r1' + a2*b2*r2';

    % compute S, z, sig, kappa, del
    S = B + B';
    sig = trace(B);
    del = det(S);

    z = [ B(2,3) - B(3,2);
            B(3,1) - B(1,3);
            B(1,2) - B(2,1) ];
    kappa = S(2,2)*S(3,3) - S(3,2)*S(2,3) + S(1,1)*S(3,3) - S(3,1)*S(1,3) + S(1,1)*S(2,2) - S(2,1)*S(1,2);

    % Compute max eigenvalue
    cRB = dot(r1,r2)*dot(b1,b2) + norm(cross(r1,r2))*norm(cross(b1,b2));
    lambda_max = sqrt(a1^2 + 2*a1*a2*cRB + a2^2);

    % Compute quaternion
    alpha = lambda_max^2 - sig^2 + kappa;
    beta = lambda_max - sig;
    gamma = (lambda_max + sig)*alpha - del;
    
    q_vec = (alpha*eye(3) + beta*S + S^2) * z;
    prefix = 1/sqrt(gamma^2 + norm(q_vec)^2);
    q_hat(:,Dataidx) = prefix * [q_vec; gamma];

    % Transform the q_hat from JPL style to Hamilton style, which is the style
    % used in MATLAB

    quat_est = quaternion(q_hat(4,Dataidx), q_hat(1,Dataidx), q_hat(2,Dataidx), q_hat(3,Dataidx));
    % Turn the quaternion into ypr angle.
    YPR(Dataidx,:) = eulerd(quat_est, "ZYX", "point");

%end



% INS ECEF ESKF
% This program is an error-state kalman filter for initial fine alignment
% using IMU data only.
%
% The IMU data imported is in ENU frame, but transformed to NED frame
% before applying the Kalman filter. 
% 
% The final output data is expressed in ECEF frame including DCM, velocity
% vector, position vector, accelerometer bias and gyroscope bias. These
% data is stored in "States" variable.
%
% The yaw, pitch and roll angle output is expressed in local NED frame
% (Navigation frame), in rad. Yawd, Pitchd and Rolld is the Euler angles in
% degree.
%
%Auther: C.T.Shen 20220902

%% Parameter definition

%%%% 
% Set this parameter to "true" to enable kalman filter correction
% Set to "false" to perform integration only.
ENABLE_KF_CORRECTION = true;
%%%%

dt = 1/200; % Sampling time
EarthRotation = 7.2921150e-5; % Earth rotation angular speed, in rad/sec.
e = 0.0818191908425;  % eccentricity of the ellipsoid.
R0 = 6378137;  % Earth's equatorial radius ,m .


TestLat = 24.72534; %24 + 43/60 + 29/3600;  % Latitude of TiSpace 
TestLong = 120.91023; %120 + 54/60 + 35/3600;  % Longitude of TiSpace 
TestAlt = 22; % 22m height above sea level
g = 9.7803253359*(1+0.001931853*sin(deg2rad(TestLat))^2)/sqrt(1-e^2*sin(deg2rad(TestLat))^2);  %m/s^2

g_LD = [0;0;g];  % Reference gravitational measurement in local NED frame.
% CHECK THE DIRECTION OF GRAVITY

w_IE_E = [0;0;EarthRotation]; % Earth angular velocity expressed in ECEF frame

%%%% Selcet the Start and End time of the testing data %%%%%%%%
Start_min = 40; % in minutes.
End_min = 60;
StartSample = 200*60*(Start_min)+1;
EndSample = 200*60*(End_min);  
if EndSample > size(accM_enu,2)
    EndSample = size(accM_enu,2);
end
DataSamples = EndSample - StartSample +1; % How many data samples used in the test

%%%% Sensor error model parameters from allan variance %%%%
SnG1 = (2.81e-06)^2;
SkG1 = (3.41e-07)^2;
SbG1 = 8.089e-11; 
TBG1 = 1;

SnG2 = (3.16e-06)^2;
SkG2 =  (1.58e-07)^2;
SbG2 = 9.38e-11; 
TBG2 = 1;

SnG3 = (2.7e-06)^2;
SkG3 = (1.89e-07)^2;
SbG3 = 4.22e-11;
TBG3 =1;

SnA1 = (2.71e-04)^2;
SkA1 = (2.18e-05)^2;
SbA1 = 1.75e-07;
TBA1 = 1;

SnA2 = (2.4e-04)^2;
SkA2 = (3.3e-05)^2;
SbA2 = 1.51e-07;
TBA2 = 1;

SnA3 = (2.04e-04)^2;
SkA3 = (5.57e-05)^2;
SbA3 = 1.99e-07;
TBA3 = 1;


%% Initialize
% Attitude matrix in ECEF
Yawd_init =  -130.5595; 
Pitchd_init = 0.0057;
Rolld_init = 179.9953;  % in NED
% LD2B = angle2dcm(deg2rad(Yawd_init),deg2rad(Pitchd_init),deg2rad(Rolld_init), 'ZYX');  % Local NED to Body (NED)
LD2B = ANGLE2DCM_321(deg2rad(Yawd_init),deg2rad(Pitchd_init),deg2rad(Rolld_init)); % Local NED to Body (NED)
% E2LD = dcmecef2ned(TestLat, TestLong);
E2LD = DCM_ECEF2NED(TestLat, TestLong);
E2B = E2LD*LD2B;
B2E_hat = E2B';

g_E = E2LD' * g_LD;  % Transform g from NED to ECEF.

v_E_hat = [0;0;0];  % Initial velocity
% r_E_hat = lla2ecef([TestLat TestLong TestAlt])';  % initial position
r_E_hat = POS_LLA2ECEF(TestLat, TestLong, TestAlt);  % initial position
bw = [-9.22E-07;3.85E-06;9.76E-07];  % initial gyro bias
ba = [-0.01173;0.01877;-0.00778]; % initial acc bias

% Pre-allocate memeory for data
States = zeros(21,DataSamples);
Yawd = zeros(1,DataSamples);
Pitchd = zeros(1,DataSamples);
Rolld = zeros(1,DataSamples);

% Initialization
States(:,1) = [ B2E_hat(:,1); B2E_hat(:,2); B2E_hat(:,3); v_E_hat; r_E_hat; ba; bw] ;
Err_states = [zeros(15,1)];
Yawd(1) = Yawd_init;
Pitchd(1) = Pitchd_init;
Rolld(1) = Rolld_init;

P_hat = 1e-6*eye(15);  % initial state covariance

Q = diag([ SnA1; SnA2; SnA3; SbA1; SbA2; SbA3; SnG1; SnG2; SnG3; SbG1; SbG2; SbG3 ]);
% Continuous process noise covariance

%% ESKF
for DataIdx = 2: DataSamples-1

    %%%% Load data and transform to NED %%%%
    wm_last = [rateM_enu(2, DataIdx+ StartSample -1) ; rateM_enu(1, DataIdx+ StartSample -1); -rateM_enu(3, DataIdx+ StartSample -1)];
    am_last = [ accM_enu(2,DataIdx + StartSample -1); accM_enu(1,DataIdx+ StartSample -1); -accM_enu(3, DataIdx+ StartSample -1)];

    wm_now = [rateM_enu(2, DataIdx+ StartSample ) ; rateM_enu(1, DataIdx+ StartSample ); -rateM_enu(3, DataIdx+ StartSample )];
    am_now = [ accM_enu(2,DataIdx + StartSample ); accM_enu(1,DataIdx+ StartSample ); -accM_enu(3, DataIdx+ StartSample )];

    
    %%%% Propagation RK4 %%%%
    % Acc Bias
    ba_k1 = -diag([1/TBA1; 1/TBA2; 1/TBA3 ])*ba;
    ba_k2 = -diag([1/TBA1; 1/TBA2; 1/TBA3 ])*(ba + (dt/2)*ba_k1);
    ba_k3 = -diag([1/TBA1; 1/TBA2; 1/TBA3 ])*(ba + (dt/2)*ba_k2);
    ba_k4 = -diag([1/TBA1; 1/TBA2; 1/TBA3 ])*(ba + dt*ba_k3);
    ba_pred = ba + dt*(ba_k1 + 2*ba_k2 + 2*ba_k3 + ba_k4)/6;

     % gyro Bias
    bw_k1 = -diag([1/TBG1; 1/TBG2; 1/TBG3 ])*bw;
    bw_k2 = -diag([1/TBG1; 1/TBG2; 1/TBG3 ])*(bw + (dt/2)*bw_k1);
    bw_k3 = -diag([1/TBG1; 1/TBG2; 1/TBG3 ])*(bw + (dt/2)*bw_k2);
    bw_k4 = -diag([1/TBG1; 1/TBG2; 1/TBG3 ])*(bw + dt*bw_k3);
    bw_pred = bw + dt*(bw_k1 + 2*bw_k2 + 2*bw_k3 + bw_k4)/6;
    
    %gyro and acc estimation
    a_pred_last = am_last - ba;
    w_pred_last = wm_last - bw;

    a_pred_now = am_now - ba_pred;
    w_pred_now = wm_now - bw_pred;

    % Attitude
    C_k1 = B2E_hat*skewMatrix(w_pred_last) - skewMatrix(w_IE_E) * B2E_hat;
    C_k2 = (B2E_hat + (dt/2)*C_k1)*skewMatrix((w_pred_now+w_pred_last)/2) - skewMatrix(w_IE_E)*(B2E_hat + (dt/2)*C_k1);
    C_k3 = (B2E_hat + (dt/2)*C_k2)*skewMatrix((w_pred_now+w_pred_last)/2) - skewMatrix(w_IE_E)*(B2E_hat + (dt/2)*C_k2);
    C_k4 = (B2E_hat + dt*C_k3)*skewMatrix(w_pred_now) - skewMatrix(w_IE_E)*(B2E_hat + dt*C_k3);
    B2E_pred = B2E_hat + dt*(C_k1 + 2*C_k2 + 2*C_k3 + C_k4)/6;
    % Specific force transformation
    a_pred_last_E = (B2E_hat*(a_pred_last - E2LD'*B2E_hat'*g_LD));  % substract gravity term at body frame
    a_pred_now_E = (B2E_pred*(a_pred_now - E2LD'*B2E_pred'*g_LD)) ;
    % Velocity
    % Because we don't know r_E_now, therefore the gravity term remains constant
    v_k1 = a_pred_last_E -skewMatrix(w_IE_E)*skewMatrix(w_IE_E)*r_E_hat - 2*skewMatrix(w_IE_E)*v_E_hat;
    v_k2 = (a_pred_now_E + a_pred_last_E)/2  -skewMatrix(w_IE_E)*skewMatrix(w_IE_E)*r_E_hat - 2 * skewMatrix(w_IE_E)*(v_E_hat+(dt/2)*v_k1);
    v_k3 = (a_pred_now_E + a_pred_last_E)/2 -skewMatrix(w_IE_E)*skewMatrix(w_IE_E)*r_E_hat - 2 * skewMatrix(w_IE_E)*(v_E_hat+(dt/2)*v_k2);
    v_k4 = a_pred_now_E -skewMatrix(w_IE_E)*skewMatrix(w_IE_E)*r_E_hat - 2 * skewMatrix(w_IE_E)*(v_E_hat+dt*v_k3);
    v_E_pred = v_E_hat + dt*(v_k1 + 2*v_k2 + 2*v_k3 + v_k4)/6;
    % Position
    r_k1 = v_E_hat;
    r_k2 = (v_E_pred + v_E_hat)/2;
    r_k3 = (v_E_pred + v_E_hat)/2;
    r_k4 = v_E_pred;
    r_E_pred = r_E_hat + dt*(r_k1 + 2*r_k2 + 2*r_k3 + r_k4)/6;
    
    
   %%%% Error state model %%%%
%     LLA_hat = ecef2lla(r_E_hat');
    LLA_hat = POS_ECEF2LLA(r_E_hat);
    Lat_hat = LLA_hat(1);
    g0 = 9.7803253359*(1+0.001931853*sin(deg2rad(Lat_hat))^2)/sqrt(1-e^2*sin(deg2rad(Lat_hat))^2);  %m/s^2
    RE = R0/sqrt(1-e^2*sin(deg2rad(Lat_hat))^2);
    r_es_E = RE*sqrt(1-e*(2-e)*sin(deg2rad(Lat_hat))^2);
    
    % Continuous-time error state state space model
    F11 = -skewMatrix(w_IE_E);
    F15 = B2E_hat;
    F21 = -skewMatrix(B2E_hat*(a_pred_last - E2LD'*B2E_hat'*g_LD));
    F22 = -2*skewMatrix(w_IE_E);
    F23 = (2*g0/ r_es_E)*(r_E_hat/norm(r_E_hat)^2)*r_E_hat';
    F24 = B2E_hat;
    F44 = -diag([1/TBA1; 1/TBA2; 1/TBA3 ]);
    F55 = -diag([1/TBG1; 1/TBG2; 1/TBG3 ]);

    G13 = B2E_hat;
    G21 = B2E_hat;

    F_INS_E = [F11 zeros(3,9) F15;
                       F21 F22 F23 F24 zeros(3,3);
                       zeros(3,3) eye(3) zeros(3,9);
                       zeros(3,9) F44 zeros(3,3);
                       zeros(3,12) F55 ];

    G_INS_E = [ zeros(3,6) G13 zeros(3,3);
                         G21 zeros(3,9);
                         zeros(3,12);
                         zeros(3,3) eye(3) zeros(3,6);
                         zeros(3,9) eye(3) ];

    % Discrete Model
    FG = [ F_INS_E G_INS_E;
                zeros(12,15) zeros(12,12) ];
    FGd = expm(FG);
    Ad = FGd(1:15,1:15);
    Bd = FGd(1:15, 16:27);

    % compute discrete process noise covariance Qd
    Aq(1:15,1:15) = -F_INS_E;
    Aq(1:15, 16:30) = G_INS_E * Q *G_INS_E';
    Aq(16:30,16:30) = F_INS_E';
    Aq = dt*Aq;
    Aqd = expm(Aq);
    Qd = Ad * Aqd(1:15, 16:30);

    Err_pred = Ad*Err_states ;

    % covariance propagation
    P_pred = Ad*P_hat*Ad' + Qd;

    %%%% Update %%%%
    
    % zero displacement
    del_pos = States(13:15,1) - r_E_pred; %hat

    % ZVU

    % Earth angularV alignment
    del_w_E = w_IE_E - B2E_pred*(w_pred_now);
    dwEdea = skewMatrix(B2E_pred*(w_pred_now));
    dwEdbw = skewMatrix(Err_pred(1:3)) * (B2E_pred);

    
    rr(:,DataIdx) = [-v_E_pred; del_pos; del_w_E];  % measurement residual

    R = diag([1e-10*ones(6,1);  1e-9*ones(3,1)]);  % R matrix


    H = [zeros(3,3) -eye(3) zeros(3,9);
            zeros(3,6) -eye(3) zeros(3,6);
            dwEdea zeros(3,9) dwEdbw ];


    K = P_pred*H' /(H*P_pred*H' + R);
    del_x = K*rr(:,DataIdx);

    Err_states =  Err_pred + del_x;

    if ENABLE_KF_CORRECTION == true
        % update B2E with error angle
        B2E_hat = (eye(3) - skewMatrix(Err_states(1:3)))*B2E_pred;
        % update velocity
        v_E_hat = v_E_pred - Err_states(4:6);
        %update position
        r_E_hat = r_E_pred - Err_states(7:9);
        %update bias
        ba =  ba_pred - Err_states(10:12);
        bw = bw_pred - Err_states(13:15);
    else
            B2E_hat = B2E_pred;
            v_E_hat = v_E_pred;
            r_E_hat = r_E_pred;
            ba = ba_pred;
            bw = bw_pred;
    end


    %update state covariance matrix 
    P_hat = (eye(15) - K*H)*P_pred*(eye(15)-K*H)' + K*R*K';

    % logging
    States(:,DataIdx) = [ B2E_hat(:,1); B2E_hat(:,2); B2E_hat(:,3); v_E_hat; r_E_hat; ba; bw] ;
    
    % Reset error state
    Err_states(1:15) = zeros(15,1);
    % clauclate new E2LD
%     LLA_hat= ecef2lla(r_E_hat');
    LLA_hat = POS_ECEF2LLA(r_E_hat);
    Lat_hat = LLA_hat(1);
    Long_hat = LLA_hat(2);
%     E2LD = dcmecef2ned(Lat_hat, Long_hat);
    E2LD = DCM_ECEF2NED(Lat_hat, Long_hat);
    LD2B_hat = E2LD'*B2E_hat';
    
    % Output Euler angle on navigation frome (NED)
%     [AA1, AA2, AA3] = dcm2angle(LD2B_hat, 'ZYX');
    [AA1, AA2, AA3] = DCM2ANGLE_321(LD2B_hat);
    Yawd(DataIdx) = rad2deg(AA1);
    Pitchd(DataIdx) = rad2deg(AA2);
    Rolld(DataIdx) = rad2deg(AA3);
    if Rolld(DataIdx) < -170
        Rolld(DataIdx) = 360 + Rolld(DataIdx);
    end

    % show progress
    if mod(DataIdx,2000) == 0
        DataIdx
    end
    

end


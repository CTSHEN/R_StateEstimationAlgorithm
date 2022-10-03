% Transform Result to ENU frame
TestLat = 24.72534; %24 + 43/60 + 29/3600;  % Latitude of TiSpace 
TestLong = 120.91023; %120 + 54/60 + 35/3600;  % Longitude of TiSpace 

E2LD = dcmecef2ned(TestLat, TestLong);
LD2LU = [0 1 0
                 1 0 0
                 0 0 -1];   % DCM from NED to ENU.

for idx = 1:2000
    B2E = [States(1:3,idx) States(4:6,idx) States(7:9,idx)];

    LU2B = E2LD' * B2E' * LD2LU';
    [RollLU(idx), PitchLU(idx), YawLU(idx)] = dcm2angle(LU2B, 'XYZ');
    RollLU(idx) = rad2deg(RollLU(idx));
    PitchLU(idx) = rad2deg(PitchLU(idx));
    YawLU(idx) = rad2deg(YawLU(idx));

    DX(idx) = States(13,idx) - States(13,1);
    DY(idx) = States(14,idx) - States(14,1);
    DZ(idx) = States(15,idx) - States(15,1);

    VX(idx) = States(10,idx);
    VY(idx) = States(11,idx);
    VZ(idx) = States(12,idx);

    BAX(idx) = States(16,idx);
    BAY(idx) = States(17,idx);
    BAZ(idx) = States(18,idx);

    BGX(idx) = States(19,idx);
    BGY(idx) = States(20,idx);
    BGZ(idx) = States(21,idx);

end

%% Plot Result
% Attitude
figure
subplot(3,1,1);
plot(RollLU);
title('Roll angle');

subplot(3,1,2);
plot(PitchLU);
title('Pitch angle');

subplot(3,1,3);
plot(YawLU);
title('Yaw angle');

% Displacement
figure
subplot(3,1,1);
plot(DX);
title('Displacement X (ECEF)');

subplot(3,1,2);
plot(DY);
title('Displacement Y (ECEF)');

subplot(3,1,3);
plot(DZ);
title('Displacement Z (ECEF)');

% Velocity
figure
subplot(3,1,1);
plot(VX);
title('Velocity X (ECEF)');

subplot(3,1,2);
plot(VY);
title('Velocity Y (ECEF)');

subplot(3,1,3);
plot(VZ);
title('Velocity Z (ECEF)');

% acc bias
figure
subplot(3,1,1);
plot(BAX);
title('AccX bias');

subplot(3,1,2);
plot(BAY);
title('AccY bias');

subplot(3,1,3);
plot(BAZ);
title('AccZ bias');

%gyro bias
figure
subplot(3,1,1);
plot(BGX);
title('GyroX bias');

subplot(3,1,2);
plot(BGY);
title('GyroY bias');

subplot(3,1,3);
plot(BGZ);
title('GyroZ bias');
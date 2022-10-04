% Plot result from start to DataIdx

%create time vector
time = zeros(1, DataIdx);
for t_idx = 1: DataIdx
    time(1, t_idx) = (1/200) * (t_idx - 1);
end

%% Choose data to be ploted
ATTITUDE = true;
VELOCITY = true;
DISPLACEMENT = true;
ACC_BIAS = true;
GYRO_BIAS = true;

%% Plot the result
if ATTITUDE == true
    figure
    subplot(3,1,1);
    plot(time, Yawd(1,1:DataIdx));
    title('Yaw (NED)');
    ylabel('degree');
    subplot(3,1,2);
    plot(time, Pitchd(1,1:DataIdx));
    title('Pitch (NED)');
    ylabel('degree');
    subplot(3,1,3);
    plot(time, Rolld(1,1:DataIdx));
    title('Roll (NED)');
    ylabel('degree');
    xlabel('time (seconds)');
end

if VELOCITY == true
    figure
    subplot(3,1,1);
    plot(time, States(10,1:DataIdx));
    title('VelocityX (ECEF)');
    ylabel('m/s');
    subplot(3,1,2);
    plot(time, States(11,1:DataIdx));
    title('VelocityY (ECEF)');
    ylabel('m/s');
    subplot(3,1,3);
    plot(time, States(12,1:DataIdx));
    title('VelocityZ (ECEF)');
    ylabel('m/s');
    xlabel('time (seconds)');
end

if DISPLACEMENT == true
    figure
    subplot(3,1,1);
    plot(time, States(13,1:DataIdx) - States(13,1));
    title('DisplacementX (ECEF)');
    ylabel('m');
    subplot(3,1,2);
    plot(time, States(14,1:DataIdx) - States(14,1));
    title('DisplacementY (ECEF)');
    ylabel('m');
    subplot(3,1,3);
    plot(time, States(15,1:DataIdx) - States(15,1));
    title('DisplacementY (ECEF)');
    ylabel('m');
    xlabel('time (seconds)');
end

if ACC_BIAS == true
    figure
    subplot(3,1,1);
    plot(time, States(16,1:DataIdx));
    title('ACCX Bias');
    ylabel('m/s^2');
    subplot(3,1,2);
    plot(time, States(17,1:DataIdx));
    title('ACCY Bias');
    ylabel('m/s^2');
    subplot(3,1,3);
    plot(time, States(18,1:DataIdx));
    title('ACCZ Bias');
    ylabel('m/s^2');
    xlabel('time (seconds)');
end

if GYRO_BIAS == true
    figure
    subplot(3,1,1);
    plot(time, States(19,1:DataIdx));
    title('GYROX Bias');
    ylabel('rad/s');
    subplot(3,1,2);
    plot(time, States(20,1:DataIdx));
    title('GYROY Bias');
    ylabel('rad/s');
    subplot(3,1,3);
    plot(time, States(21,1:DataIdx));
    title('GYROZ Bias');
    ylabel('rad/s');
    xlabel('time (seconds)');
end
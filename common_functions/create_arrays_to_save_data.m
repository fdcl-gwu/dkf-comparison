function [true, imu, vicon, est_att, est_pos, gps, lidar] = ...
    create_arrays_to_save_data(N, m, n)

% Simulated true data
true.x = zeros(3, N);
true.v = zeros(3, N);
true.a = zeros(3, N);
true.R = zeros(3, 3, N);
true.W = zeros(3, N);
true.b_a = zeros(1, N);

% Estimated data
est_pos.x = zeros(3, N);
est_pos.v = zeros(3, N);
est_pos.a = zeros(3, N);
est_pos.b_a = zeros(1, N);
est_pos.P = zeros(m, m, N);
est_pos.R = zeros(3, 3, N);
est_pos.W = zeros(3, N);

est_att.R = zeros(3, 3, N);
est_att.W = zeros(3, N);
est_att.P = zeros(n, n, N);

% IMU measurements
imu.R = zeros(3, 3, N);
imu.a = zeros(3, N);
imu.W = zeros(3, N);

% Vicon measurements
vicon.R = zeros(3, 3, N);
vicon.x = zeros(3, N);

% GPS measurements
gps.x = zeros(3, N);
gps.v = zeros(3, N);

% Lidar measurements
lidar.d = zeros(1, N);

end
close all;
clear;

addpath('common_functions/');  % includes hat, vee, etc
%% Define variables
tf = 20;  % final time
freq = 100;  % frequency of the system (Hz)
freq_gps = 5;  % frequency of the GPS update (Hz)

% Parameters
R_bi = eye(3);

gps_delay = 0.2;  % GPS measurement delay (this code is hard coded
    % such that this delay must be a multiplication of freq_gps)

% Variances of w_k
V_a = 1e-1*diag([0.1 0.1 0.1]).^2;  % acceleration
V_W = 1e-1*diag([0.5 0.5 0.5]).^2;  % angular velocity
V_b_a = 0.01^2;  % acclerometer z bias

Q = blkdiag(V_a, V_W, V_b_a);

% Covariances of measurement error
V_R_imu = diag([0.01 0.01 0.01]).^2;
V_x_gps = diag([0.01 0.01 0.01]).^2;
V_v_gps = diag([0.01 0.01 0.01]).^2;

V_x_lidar = 0.01^2;

% Initial covariances of x
P_x = diag([1^2, 1^2, 1^2]);  % position
P_v = diag([1^2, 1^2, 1^2]);  % velocity
P_eta = diag([0.1^2, 0.1^2, 0.1^2]);  % attitude
P_b_a = 1^2;  % accelerometer z bias

P = blkdiag(P_x, P_v, P_eta, P_b_a);
[m, ~] = size(P);

%% Calculate time related parameters
N = tf*freq + 1;
t = linspace(0,tf,N);
h = t(2) - t(1);

%% Create empty arrays to save data
% tru: simulated true data
% imu: simulated IMU measurement with noise
% vcn: simulated Vicon data with noise
% est: estimated data
% gps: simulated GPS data with noise
% ldr: simulated Lidar data with noise
[tru, imu, vcn, est, gps, ~] = create_arrays_to_save_data(N, m);

%% Initial estimates
est.x(:,1) = [0, 0, 0]';
est.v(:,1) = [0, 0, 0]';
est.R(:,:,1) = eye(3);
est.b_a(1) = 0;
est.P(:,:,1) = P;

DKF = DelayedKalmanFilter();
DKF.update_parameters(R_bi);
DKF.update_init_values(P, Q);
%% Run estimator
[tru.x(:,1), tru.v(:,1), tru.a(:,1), tru.R(:,:,1), ...
    tru.W(:,1), tru.b_a(:,1)] = true_data(t(1));
[imu.a(:,1), imu.R(:,:,1), imu.W(:,1)] = measurement_imu(t(1), ...
        R_bi, V_R_imu, V_W, V_a);

% Based on the IMU and GPS requencies, we should expect new GPS readings
% at every following loop.
k_gps = round(freq / freq_gps);
N_gps = k_gps/freq_gps;

norm_error = 0;

for k = 2:N-1-k_gps
    % Update simulated true data and measurements
    [tru.x(:,k), tru.v(:,k), tru.a(:,k), tru.R(:,:,k), ...
        tru.W(:,k), tru.b_a(:,1)] = true_data(t(k));
    [imu.a(:,k), imu.R(:,:,k), imu.W(:,k)] = measurement_imu(t(k), ...
        R_bi, V_R_imu, V_W, V_a);
    
    gps_new_data = gps_measurement_check(k, k_gps);
    
    if gps_new_data
        % No delay
        % [gps.x(:,k), gps.v(:,k)] = measurement_gps(t(k), ...
        %     V_x_gps, V_v_gps);
        
        % GPS measurements are delayed by k_gps steps.
        [gps.x(:,k), gps.v(:,k)] = measurement_gps(t(k - k_gps), ...
            V_x_gps, V_v_gps);
    else
        gps.x(:,k) = gps.x(:,k-1);
        gps.v(:,k) = gps.v(:,k-1);
    end
    
    % Run prediction
    DKF.prediction(h, imu.a(:,k), imu.W(:,k));
    
    % GPS correction
    if gps_new_data
        % No delay
        % DKF.correction_gps(gps.x(:,k), gps.v(:,k), V_x_gps, V_v_gps);
        
        % Delayed GPS data
        k_list = linspace(k-k_gps, k, k_gps/N_gps);
        % k_list = k-k_gps:k; % re-calculate
        
        DKF.start_delayed_gps_correction(est.x(:,k-k_gps), ...
            est.v(:,k-k_gps), est.a(:,k-k_gps), est.b_a(k-k_gps), ...
            est.R(:,:,k-k_gps), est.W(:,k-k_gps), ...
            est.P(:,:,k-k_gps), ...
            imu.a(:,k-k_gps), imu.W(:,k-k_gps));
        DKF.correction_delayed_gps(gps.x(:,k), gps.v(:,k), ...
            V_x_gps, V_v_gps);
        
        for i = 1:length(k_list)-1
            k1 = k_list(i);
            k2 = k_list(i+1);
            h_s = t(k2) - t(k1);
            DKF.forward_prediction(h_s, imu.a(:,k2), imu.W(:,k2));
        end
        DKF.end_delayed_gps_correction();
    end
    
    % IMU correction
    DKF.correction_imu(imu.R(:,:,k), V_R_imu);
    
    [est.x(:,k), est.v(:,k), est.a(:,k), est.b_a(k), ...
        est.R(:,:,k), est.W(:,k),  est.P(:,:,k)] ...
        = DKF.output_states();
    
    norm_error = norm_error + norm(est.x(:,k) - tru.x(:,k));
end

%% Plots

plot_3x1_3e(t, tru.x, gps.x, est.x,  "$x$ (m)");  
plot_3x1_3e(t, tru.v, gps.v, est.v, "$v$ (m/s)");

%% Create empty arrays to save data, reduce clutter int the front
function [true, imu, vicon, est, gps, lidar] = ...
    create_arrays_to_save_data(N, m)

% Simulated true data
true.x = zeros(3, N);
true.v = zeros(3, N);
true.a = zeros(3, N);
true.R = zeros(3, 3, N);
true.W = zeros(3, N);

% Estimated data
est.x = zeros(3, N);
est.v = zeros(3, N);
% est.a = zeros(3, N);
est.R = zeros(3, 3, N);
est.W = zeros(3, N);
est.b_a = zeros(1, N);
est.P = zeros(m, m, N);

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
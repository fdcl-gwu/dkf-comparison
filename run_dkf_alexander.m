close all;
clear;

addpath('common_functions/');  % includes hat, vee, etc

%% Define variables
tf = 20;  % final time
freq = 100;  % frequency of the system (Hz)
freq_gps = 5;  % frequency of the GPS update (Hz)
gps_delay = 0.2;  % GPS measurement delay (this code is hard coded
    % such that this delay must be a multiplication of freq_gps)
    
N = tf * freq + 1;
t = linspace(0, tf, N);
h = t(2) - t(1);

% Parameters
R_bi = eye(3);
g = 9.81;
%% Covariances
% Variances of w_k
V_a = 1e-1*diag([0.1 0.1 0.1]).^2;  % acceleration
V_W = 1e-1*diag([0.5 0.5 0.5]).^2;  % angular velocity
V_b_a = 0.05^2;  % acclerometer z bias
V_R_imu = diag([0.01 0.01 0.01]).^2;

V_x_gps = 0.01^2 * eye(3);
V_v_gps = 0.01^2 * eye(3);
V_a_imu = 0.01^2 * eye(3);

% Initial covariances of x
P_x = diag([1^2, 1^2, 1^2]);  % position
P_v = diag([1^2, 1^2, 1^2]);  % velocity
P_eta = diag([0.1^2, 0.1^2, 0.1^2]);  % attitude
P_b_a = 0.1^2;  % accelerometer z bias

P_att0 = P_eta;
Q_att = V_a;

P_pos0 = blkdiag(P_x, P_v, P_b_a);
Q_pos = blkdiag(V_x_gps, V_v_gps, V_b_a);

%% Init filter
DKF = DelayedKalmanFilterAlexandar();
DKF.update_parameters(R_bi);
DKF.update_init_values(P_att0, P_pos0, Q_att, Q_pos)

%% Create empty arrays to save data
% tru: simulated true data
% imu: simulated IMU measurement with noise
% vcn: simulated Vicon data with noise
% est: estimated data
% gps: simulated GPS data with noise
% ldr: simulated Lidar data with noise
[tru, imu, vcn, est_att, est_pos, gps, ldr] = ...
    create_arrays_to_save_data(N, 7, 3);

%% Run estimator
[tru.x(:,1), tru.v(:,1), tru.a(:,1), tru.R(:,:,1), ...
    tru.W(:,1), tru.b_a(1)] = true_data(t(1));
[imu.a(:,1), imu.R(:,:,1), imu.W(:,1)] = measurement_imu(t(1), ...
    R_bi, V_R_imu, V_W, V_a);

% Based on the IMU and GPS requencies, we should expect new GPS readings
% at every following loop.
k_gps = round(freq/freq_gps);
N_gps = k_gps/freq_gps;

norm_error = 0;

t_prev = 0;
for k = 2:N-1-k_gps
    % Update simulated true data and measurements
    [tru.x(:,k), tru.v(:,k), tru.a(:,k), tru.R(:,:,k), ...
        tru.W(:,k), tru.b_a(k)] = true_data(t(k));
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
        
    % Run predictions
    DKF.prediction_position(h, imu.a(:,k));
    DKF.prediction_attitude(h, imu.W(:,k));
    
    % GPS correction
    if gps_new_data
        % No delay
        % DKF.correction_gps(gps.x(:,k), gps.v(:,k), V_x_gps, V_v_gps);

        k_list = linspace(k-k_gps,k,k_gps/N_gps);
        
        % Delayed GPS data
        DKF.start_delayed_gps_correction(est_pos.x(:,k-k_gps), ...
            est_pos.v(:,k-k_gps), est_pos.b_a(k-k_gps));
        DKF.correction_delayed_gps(gps.x(:,k), gps.v(:,k), ...
            V_x_gps, V_v_gps);
        DKF.end_delayed_gps_correction();
        
        % Correct P as if there is no delay, i.e. when the measurement
        % was supposed to arrive, assuming H and R are known.
        % This is done here because the GPS delay and the
        % measurement gap are the same.
        DKF.postion_correction_P(V_x_gps, V_v_gps);
    end
    
    % IMU correction
    DKF.correction_imu(imu.R(:,:,k), V_R_imu);
    
    [est_pos.x(:,k), est_pos.v(:,k), est_pos.a(:,k), est_pos.b_a(k), ...
        est_att.R(:,:,k), est_att.W(:,k), ...
        est_att.P(:,:,k), est_pos.P(:,:,k)] = DKF.output_states();
    
    norm_error = norm_error + norm(est_pos.x(:,k) - tru.x(:,k));
end

%% Plots

plot_3x1_3e(t, tru.x, gps.x, est_pos.x,  "$x$ (m)");  
plot_3x1_3e(t, tru.v, gps.v, est_pos.v, "$v$ (m/s)");
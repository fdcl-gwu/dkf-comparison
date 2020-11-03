close all;
clear;

addpath('common_functions/');  % includes hat, vee, etc

%% Define variables
tf = 20;  % final time
freq = 100;  % frequency of the system (Hz)
freq_gps = 5;  % frequency of the GPS update (Hz)
gps_delay = 0.2;  % GPS measurement delay (this code is hard coded
    % such that this delay must be a multiplication of freq_gps)

is_delayed = true;  % flag to change delayed/non-delayed mode for GPS
    
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
V_b_a = 0.001^2;  % acclerometer z bias
V_R_imu = diag([1, 1, 1])*0.01^2;

V_x_gps = 0.001^2 * eye(3);
V_v_gps = 0.001^2 * eye(3);
V_a_imu = 0.001^2 * eye(3);

% Initial covariances of x
P_x = diag([1^2, 1^2, 1^2]);  % position
P_v = diag([1^2, 1^2, 1^2]);  % velocity
P_eta = diag([0.1^2, 0.1^2, 0.1^2]);  % attitude
P_b_a = 0.01^2;  % accelerometer z bias

Q_att = V_a;

P = blkdiag(P_x, P_v, P_eta, P_b_a);
Q = blkdiag(V_x_gps, V_v_gps, V_R_imu, V_b_a);

%% Init filter
DKF = DelayedKalmanFilterMerwe();
DKF.update_parameters(R_bi);
DKF.update_init_values(P, Q);

%% Create empty arrays to save data
% tru: simulated true data
% imu: simulated IMU measurement with noise
% vcn: simulated Vicon data with noise
% est: estimated data
% gps: simulated GPS data with noise
% ldr: simulated Lidar data with noise
[tru, imu, vcn, ~, est, gps, ldr] = ...
    create_arrays_to_save_data(N, 7, 3);
est.P = zeros(10, 10, N);

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
    ypr = rotm2eul(imu.R(:,:,k))';
  
    gps_new_data = gps_measurement_check(k, k_gps);

    if gps_new_data
        if is_delayed
            % GPS measurements are delayed by k_gps steps.
            [gps.x(:,k), gps.v(:,k)] = measurement_gps(t(k - k_gps), ...
                V_x_gps, V_v_gps);
        else
            % No delay
            [gps.x(:,k), gps.v(:,k)] = measurement_gps(t(k), ...
                V_x_gps, V_v_gps);
        end
    else
        gps.x(:,k) = gps.x(:,k-1);
        gps.v(:,k) = gps.v(:,k-1);

    end
        
    % Run predictions
    if is_delayed
        DKF.delayed_prediction(h, imu.a(:,k), imu.W(:,k));
    else
        DKF.prediction(h, imu.a(:,k), imu.W(:,k));
    end
    
    % Run corrections
    if is_delayed
        DKF.delayed_correction_imu(ypr, V_R_imu);
    else
        DKF.correction_imu(ypr, V_R_imu);
    end
    
    if gps_new_data 
        % GPS correction
        if is_delayed
            DKF.delayed_correction_gps(gps.x(:,k), gps.v(:,k), V_x_gps, ...
                V_v_gps);
            DKF.update_augmented_state();
        else
            DKF.correction_gps(gps.x(:,k), gps.v(:,k), V_x_gps, V_v_gps);
        end
    end
    
    [est.x(:,k), est.v(:,k), est.R(:,:,k), est.P(:,:,k)] ...
        = DKF.output_states();
    
    norm_error = norm_error + norm(est.x(:,k) - tru.x(:,k));
end
%% Plots
plot_3x1_3e(t, tru.x, gps.x, est.x,  "$x$ (m)");  
plot_3x1_3e(t, tru.v, gps.v, est.v, "$v$ (m/s)");
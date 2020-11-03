function [a_i, R_ni, W_i] = measurement_imu(t, R_bi, ...
    V_I_R, Q_W, Q_a)

e3 = [0, 0, 1]';

% Get true data, add noise to them and simulate as sensor data.
[~, ~, a, R_fb, W_b, b_a] = true_data(t);

% Acceleration
% Accelerometers measure gravity as positive accleration upwards earth.
w_a = mvnrnd(zeros(3,1), Q_a)';
a_i = R_bi'*R_fb'*(a - b_a*e3) + w_a;
% a_i = R_bi'*R_fb'*(a - b_a) + w_a;

% Angular velocity
w_W = mvnrnd(zeros(3,1), Q_W)';
W_i = R_bi'*W_b + w_W;

% Attitude
zeta = mvnrnd(zeros(3,1), V_I_R)';
R_ni = R_fb*R_bi*expm_SO3(zeta);

end
function [x_f, v_f] = measurement_gps(t, V_x_gps, V_v_gps)

% Get true data, add noise to them and simulate as sensor data.
[x_g, v_g, ~, ~, ~] = true_data(t);

% Position
zeta_x = mvnrnd(zeros(3, 1), V_x_gps)';
x_f = x_g + zeta_x;

% Velocity
zeta_v = mvnrnd(zeros(3, 1), V_v_gps)';
v_f = v_g + zeta_v;

end
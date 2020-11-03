function [x, v, a, R, W, b_a] = true_data(t)

w1 = 0.2;
w2 = 0.1;
x = [1.2*sin(w1*pi*t), 4.2*cos(w2*pi*t), -0.5*t]';
v = [1.2*w1*pi*cos(w1*pi*t), -4.2*w2*pi*sin(w2*pi*t), -0.5]';
a = [-1.2*(w1*pi)^2*sin(w1*pi*t), -4.2*(w2*pi)^2*cos(w2*pi*t), 0]';
b_a = 1.5;
% b_a = 1.5*[1, 1, 1]';


% R = [cos(t), -cos(t)*sin(t), sin(t)^2;
%     cos(t)*sin(t), cos(t)^3 - sin(t)^2, -cos(t)*sin(t) - cos(t)^2*sin(t);
%     sin(t)^2, cos(t)*sin(t) + cos(t)^2*sin(t), cos(t)^2 - cos(t)*sin(t)^2];
% W = [cos(t) + 1, sin(t) - sin(2*t)/2, cos(t) - cos(t)^2 + 1]';

w = 0.1;
wt = w * t;
R = [cos(wt), -cos(wt)*sin(wt), sin(wt)^2;
    cos(wt)*sin(wt), cos(wt)^3 - sin(wt)^2, -cos(wt)*sin(wt) - cos(wt)^2*sin(wt);
    sin(wt)^2, cos(wt)*sin(wt) + cos(wt)^2*sin(wt), cos(wt)^2 - cos(wt)*sin(wt)^2];
W = w*[cos(wt) + 1, sin(wt) - sin(2*wt)/2, cos(wt) - cos(wt)^2 + 1]';

% Limit roll and pitch such that lidar sensor is always oriented towards
% the floor, comment if Lidar is not used.
% t = t / 50;
% th1 = t / 180 * pi;
% R1 = [1, 0, 0;
%     0, cos(th1), -sin(th1);
%     0, sin(th1), cos(th1)];
% 
% th2 = t / 180 * pi;
% R2 = [cos(th2), 0, sin(th2);
%     0, 1, 0;
%     -sin(th2), 0, cos(th2)];
% R = R2 * R1;
% W = [pi / 180, pi / 180, 0]';


% TEMP
% x = [1*t, 0, 1*t]';
% v = [1, 0, 1]';
% x = [0, 0, 0]';
% v = [0, 0, 0]';
% a = [0, 0, 0]';
b_a = 0;
% 
% R = eye(3);
% W = [0, 0, 0]';

end
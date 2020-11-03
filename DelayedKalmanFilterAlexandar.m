classdef DelayedKalmanFilterAlexandar < handle
    %%
    properties (SetAccess = private)
        % States
        x
        v
        a
        b_a
        R
        W
        
        P_att
        Q_att
        
        P_pos
        Q_pos
        
        % Delayed states
        x_s
        v_s
        b_a_s
        P_s
        M
        
        W_pre
        a_imu_pre
        R_pre;
        
        R_bi
        R_fv
        R_nv
        
        new_cycle
    end
    %%
    properties (Constant)
        g = 9.81;
        e3 = [0, 0, 1]';
        ge3 = 9.81*[0, 0, 1]';
        I3 = eye(3);
        Z3 = zeros(3);
    end
    %%
    methods
        %%
        function E = DelayedKalmanFilterAlexandar()
            E.x = [0, 0, 0]';
            E.v = [0, 0, 0]';
            E.a = [0, 0, 0]';
            E.b_a = 0;
            E.R = eye(3);
            E.W = [0, 0, 0]';
            
            E.P_att = zeros(3);
            E.Q_att = zeros(3);
            
            E.P_pos = zeros(7);
            E.Q_pos = zeros(7);
            
            E.x_s = [0, 0, 0]';
            E.v_s = [0, 0, 0]';
            E.b_a_s = 0;
            E.M = zeros(7);
            
            E.R_bi = eye(3);
            
            E.W_pre = [0, 0, 0]';
            E.a_imu_pre = [0, 0, 0]';
            E.R_pre = eye(3);
            
            E.new_cycle = true;
        end
        %% Attitude Prediction
        function prediction_attitude(E, h, W_imu)
            E.W = E.R_bi*W_imu;
            E.R = E.R*expm_SO3(h/2*(E.W + E.W_pre));
            
            % Calculate A(t_{k-1})
            A = -hat(E.W_pre);
            
            % Calculate F(t_{k-1})
            F = E.R_bi;
            
            % Calculate \Psi using A(t)
            psi = E.I3 + h/2*A ;
            
            A = E.I3 + h*A*psi;
            F = h*psi*F;
            
            E.P_att = A*E.P_att*A' + F*E.Q_att*F';
            
            E.W_pre = E.W;
        end
        %% Position Prediction
        function prediction_position(E, h, a_imu)
            A = [E.I3, h*E.I3, h^2*E.e3; 
                E.Z3, E.I3, h*E.e3;
                zeros(1, 3), zeros(1, 3), 1];
            B = [h^2*E.I3; h*E.I3; zeros(1, 3)];
            u = E.R_pre*E.R_bi*E.a_imu_pre;
            
            x_k_1 = [E.x; E.v; E.b_a];
            x_k = A*x_k_1 + B*u;
            
            E.x = x_k(1:3);
            E.v = x_k(4:6);
            E.b_a = x_k(7);
            
            E.P_pos = A*E.P_pos*A' + E.Q_pos;
            
            E.a_imu_pre = a_imu;
            E.R_pre = E.R;
        end
        %% Correction from IMU measurements
        function correction_imu(E, R_imu, V_R_imu)
            imu_R = E.R'*R_imu*E.R_bi';
            del_z = 0.5*vee(imu_R - imu_R');
            
            H = E.I3;
            G = E.R_bi;
            V = V_R_imu;
            
            S = H*E.P_att*H' + G*V*G';
            K = E.P_att*H'*inv(S);
            X = K*del_z;
            
            eta = X;
            E.R = E.R*expm_SO3(eta);
            
            I_KH = E.I3 - K * H;
            E.P_att = I_KH*E.P_att*I_KH' + K*G*V*G'*K';
        end
        %% Position Correction P
        function postion_correction_P(E, V_x_gps, V_v_gps)
            H = [E.I3, E.Z3, zeros(3, 1);
                E.Z3, E.I3, zeros(3, 1)];
            V = [V_x_gps, E.Z3;
                E.Z3, V_v_gps];
            
            S = H*E.P_pos*H' + V;
            K = E.P_pos*H'*inv(S);
            
            I_KH = eye(7) - K*H;
            E.P_pos = I_KH*E.P_pos*I_KH' + K*V*K';
        end
        %% Correction from delayed GPS measurements
        function correction_delayed_gps(E, x_gps, v_gps, V_x_gps, V_v_gps)
            H = [E.I3, E.Z3, zeros(3, 1); 
                E.Z3, E.I3, zeros(3, 1)];
            V = [V_x_gps, zeros(3);
                V_v_gps, zeros(3)];
            
            S = H*E.P_pos*H' + V;
            K = E.P_pos*H'*inv(S);
            
            x_k_s = [E.x_s; E.v_s; E.b_a_s];
            x_k = [E.x; E.v; E.b_a];
            
            z = [x_gps; v_gps];
            del_x_hat = E.M*K*(z - H*x_k_s);
            
            x_k = x_k + del_x_hat;
            
            E.x = x_k(1:3);
            E.v = x_k(4:6);
            E.b_a = x_k(7);
        end
        %% Start delayed GPS correction
        % This must be called before running correction_delayed_gps
        function start_delayed_gps_correction(E, x, v, b_a)
            E.x_s = x;
            E.v_s = v;
            E.b_a_s = b_a;
        end
        %% End delayed GPS correction
        % This must be called after running forward_prediction
        function end_delayed_gps_correction(E)
            E.M = eye(7);
            E.new_cycle = true;
        end
        %% Update parameters
        function update_parameters(E, R_bi)
            E.R_bi = R_bi;
        end
        %% Update initial values
        function update_init_values(E, P_att, P_pos, Q_att, Q_pos)
            E.P_att = P_att;
            E.P_pos = P_pos;
            E.Q_att = Q_att;
            E.Q_pos = Q_pos;
        end
        %% Get States
        function [x, v, a, b_a, R, W, P_att, P_pos] = output_states(E)
            x = E.x;
            v = E.v;
            a = E.a;
            b_a = E.b_a;
            R = E.R;
            W = E.W;
            P_att = E.P_att;
            P_pos = E.P_pos;
        end
    end
end
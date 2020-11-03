classdef DelayedKalmanFilter < handle
    %%
    properties (SetAccess = private)
        % States
        x
        v
        a
        b_a
        R
        W
        P
        Q
        
        % Delayed states
        x_s
        v_s
        b_a_s
        a_s
        R_s
        W_s
        P_s
        
        a_s_pre
        b_a_s_pre
        R_s_pre
        W_s_pre
        a_imu_s_pre
        W_imu_s_pre
        
        W_pre
        a_imu_pre
        R_pre
        b_a_pre
        
        R_bi
        R_fv
        R_nv
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
        function E = DelayedKalmanFilter()
            E.x = [0, 0, 0]';
            E.v = [0, 0, 0]';
            E.a = [0, 0, 0]';
            E.b_a = 0;
            E.R = E.I3;
            E.W = [0, 0, 0]';
            E.P = zeros(10);
            E.Q = eye(9);
            
            E.x_s = [0, 0, 0]';
            E.v_s = [0, 0, 0]';
            E.a_s = [0, 0, 0]';
            E.b_a_s = 0;
            E.R_s = eye(3);
            E.W_s = [0, 0, 0]';
            E.P_s = zeros(10);
            
            E.R_bi = E.I3;
            
            E.W_pre = [0, 0, 0]';
            E.a_imu_pre = [0, 0, 0]';
            E.R_pre = E.I3;
            E.b_a_pre = 0;
            
            E.a_s_pre = [0, 0, 0]';
            E.b_a_s_pre = 0;
            E.R_s_pre = E.I3;
            E.W_s_pre = [0, 0, 0]';
            E.a_imu_s_pre = [0, 0, 0]';
            E.W_imu_s_pre = [0, 0, 0]';
        end
        %% Prediction
        function prediction(E, h, a_imu, W_imu)
            E.R_pre = E.R;
            E.W_pre = E.W;
            E.b_a_pre = E.b_a;
            
            E.W = E.R_bi*W_imu;
            E.R = E.R*expm_SO3(h/2*(E.W + E.W_pre));
            
            % This assumes IMU provide acceleration without g
            E.a = E.R*E.R_bi*a_imu + E.b_a*E.e3;
            a_pre = E.R_pre*E.R_bi*E.a_imu_pre + E.b_a_pre*E.e3;
               
            E.x = E.x + h*E.v + h^2/2*a_pre;
            E.v = E.v + h/2*(E.a + a_pre);
            
            % Calculate A(t_{k-1})
            A = zeros(10);
            A(1:3,4:6) = E.I3;
            A(4:6,7:9) = -E.R_pre*hat(E.R_bi*E.a_imu_pre);
            A(4:6,10) = E.e3;
            A(7:9,7:9) = -hat(E.R_bi*W_imu);
            
            % Calculate F(t_{k-1})
            F = zeros(10, 7);
            F(4:6,1:3) = E.R_pre*E.R_bi;
            F(7:9,4:6) = E.R_bi;
            F(10,7) = 1;
            
            % Calculate \Psi using A(t)
            psi = eye(10) + h/2*A ;
            
            A = eye(10) + h*A*psi;
            F = h*psi*F;
            
            E.P = A*E.P*A' + F*E.Q*F';
            
            E.a_imu_pre = a_imu;
        end
        %% Correction from IMU measurements
        function correction_imu(E, R_imu, V_R_imu)
            imu_R = E.R'*R_imu*E.R_bi';
            del_z = 0.5*vee(imu_R - imu_R');
            
            H = [E.Z3, E.Z3, E.I3, zeros(3, 1)];
            G = E.R_bi;
            V = V_R_imu;
            
            S = H*E.P*H' + G*V*G';
            K = E.P*H'*inv(S);
            
            X = K*del_z;
            
            eta = X(7:9);
            E.R = E.R*expm_SO3(eta);
            
            I_KH = eye(10) - K*H;
            E.P = I_KH*E.P*I_KH' + K*G*V*G'*K';
        end
        %% Correction from GPS measurements
        function correction_gps(E, x_gps, v_gps, V_x_gps, V_v_gps)   
            del_z = [x_gps - E.x;
                v_gps - E.v];
            
            H = [E.I3, E.Z3, E.Z3, zeros(3, 1);
                E.Z3, E.I3, E.Z3, zeros(3, 1);];
            G = [E.I3, E.Z3;
                E.Z3, E.I3];
            V = [V_x_gps, E.Z3;
                E.Z3, V_v_gps];
            
            S = H*P*H' + G*V*G';
            K = P*H'*inv(S);
            X = K*del_z;
            
            dx = X(1:3);
            dv = X(4:6);
            db_a = X(10);
            
            E.x = E.x + dx;
            E.v = E.v + dv;
            E.b_a = E.b_a + db_a;

            I_KH = eye(10) - K*H;
            E.P = I_KH*E.P*I_KH' + K*G*V*G'*K';
        end
        %% Correction from delayed GPS measurements
        function correction_delayed_gps(E, x_gps, v_gps, V_x_gps, V_v_gps)
            del_z = [x_gps - E.x_s;
                v_gps - E.v_s];
            
            H = [E.I3, E.Z3, E.Z3, zeros(3, 1);
                E.Z3, E.I3, E.Z3, zeros(3, 1);];
            G = [E.I3, E.Z3;
                E.Z3, E.I3];
            V = [V_x_gps, E.Z3;
                E.Z3, V_v_gps];
            
            S = H*E.P*H' + G*V*G';
            K = E.P*H'*inv(S);
            
            X = K * del_z;
            
            dx = X(1:3);
            dv = X(4:6);
            db_a = X(10);
            
            E.x_s = E.x_s + dx;
            E.v_s = E.v_s + dv;
            E.b_a_s = E.b_a_s + db_a;
            
            I_KH = eye(10) - K*H;
            E.P_s = I_KH*E.P_s*I_KH' + K*V*K';
        end
        %% Forward propagation
        function forward_prediction(E, h, a_imu_s, W_imu_s)
            E.a_s_pre = E.a_s;
            E.b_a_s_pre = E.b_a_s;
            E.R_s_pre = E.R_s;
            E.W_s_pre = E.W_s;
            
            E.W_s = E.R_bi*W_imu_s;
            E.R_s = E.R_s*expm_SO3(h/2*(E.W_s + E.W_s_pre));
            
            E.a_s = E.R_s*E.R_bi*a_imu_s + E.b_a_s*E.e3;

            E.x_s = E.x_s + h*E.v_s + h^2/2*E.a_s_pre;
            E.v_s = E.v_s + h/2*(E.a_s + E.a_s_pre);
            
            % Calculate A(t_{k-1})
            A = zeros(10);
            A(1:3,4:6) = eye(3);
            A(4:6,7:9) = -E.R_s_pre*hat(E.R_bi*E.a_imu_s_pre);
            A(4:6,10) = E.e3;
            A(7:9,7:9) = -hat(E.R_bi*E.W_imu_s_pre);
            
            % Calculate F(t_{k-1})
            F = zeros(10, 7);
            F(4:6,1:3) = E.R_s_pre*E.R_bi;
            F(7:9,4:6) = E.R_bi;
            F(10,7) = 1;
            
            % Calculate \Psi using A(t)
            psi = eye(10) + h/2*A;
            
            A = eye(10) + h*A*psi;
            F = h *psi*F;
            
            E.P_s = A*E.P_s*A' + F*E.Q*F';
            
            E.a_imu_s_pre = a_imu_s;
            E.W_imu_s_pre = W_imu_s;
        end
        %% Start delayed GPS correction
        % This must be called before running correction_delayed_gps
        function start_delayed_gps_correction(E, x, v, a, b_a, R, W, P, ...
                a_imu, W_imu)
            E.x_s = x;
            E.v_s = v;
            E.a_s = a;
            E.R_s = R;
            E.W_s = W;
            E.b_a_s = b_a;
            E.P_s = P;
            
            E.a_imu_s_pre = a_imu;
            E.W_imu_s_pre = W_imu;
        end
        %% End delayed GPS correction
        % This must be called after running forward_prediction
        function end_delayed_gps_correction(E)
            E.x = E.x_s;
            E.v = E.v_s;
            E.a = E.a_s;
            E.b_a = E.b_a_s;
            E.P = E.P_s;
        end
        %% Update parameters
        function update_parameters(E, R_bi)
            E.R_bi = R_bi;
        end
        %% Update initial values
        function update_init_values(E, P, Q)
            E.P = P;
            E.Q = Q;
        end
        %% Get States
        function [x, v, a, b_a, R, W, P] = output_states(E)
            x = E.x;
            v = E.v;
            a = E.a;
            b_a = E.b_a;
            R = E.R;
            W = E.W;
            P = E.P;
        end
    end
end
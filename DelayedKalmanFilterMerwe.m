classdef DelayedKalmanFilterMerwe < handle
    %%
    properties (SetAccess = private)
        % States
        x
        
        P
        Q
        
        n
        m
        
        kappa
        alpha
        lambda
        gamma
        
        Xi
        Wi
        
        % Delayed states
        x_s
        P_s
        Q_s
        
        Xi_s
        Wi_s
        n_s
        m_s
        
        a
        
        W_pre
        a_imu_pre
        
        R_bi
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
        function E = DelayedKalmanFilterMerwe()
            
            n = 10;
            m = 2*n + 1;
            % alpha = 1.0;
            alpha = 1.45;
            kappa = 10;
            
            E.n = n;
            E.m = m;
            E.kappa = kappa;
            E.alpha = alpha;
            
            E.lambda = alpha^2 * (n + kappa) - n;
            E.gamma = sqrt(n + E.lambda);
            
            E.x = zeros(n, 1);
            E.Xi = zeros(n, m);
            E.Wi = zeros(m, 1);
            
            E.P = zeros(n);
            
            E.n_s = 2*n;
            E.m_s = 2*E.n_s + 1;
            E.x_s = zeros(E.n_s, 1);
            E.P_s = zeros(E.n_s);
            E.Xi_s = zeros(E.n_s, E.m_s);
            E.Wi_s = zeros(E.m_s, 1);
            
            E.R_bi = eye(3);
            
            E.a = [0, 0, 0]';
            E.a_imu_pre = [0, 0, 0]';
            E.W_pre = [0, 0, 0]';
        end
        %% Prediction
        function prediction(E, h, a_imu, W_imu)
            
            E.sigma_points();
            
            Xk = zeros(E.n, E.m);
            for i = 1:E.m
                Xk(:,i) = f(E, E.Xi(:, i), h, a_imu, W_imu);
            end
            E.Xi = Xk;
            
            [E.x, E.P] = unscented_transformation(E, Xk, E.Wi, E.Q);
        end
        %% IMU correction
        function correction_imu(E, ypr, V_R_imu)
            
            % z = reshape(R, 9, 1);
            z = ypr;
            R = V_R_imu;
            
            Y = zeros(3, E.m);
            for i = 1:E.m
                Y(:,i) = E.Xi(7:9,i);
            end
            
            [y, Py] = unscented_transformation(E, Y, E.Wi, R);
            
            Pxz = zeros(E.n, 3);
            for i = 1:E.m
                Pxz = Pxz + E.Wi(i)*(E.Xi(:,i) - E.x)*(Y(:,i) - y)';
            end
            
            K = Pxz*inv(Py);
            
            E.x = E.x + K*(z - y);
            E.P = E.P - K*Py*K';
        end
        %% GPS correction
        function correction_gps(E, x_gps, v_gps, V_x_gps, V_v_gps)
            
            z = [x_gps; v_gps];
            R = blkdiag(V_x_gps, V_v_gps);
            
            Y = zeros(6, E.m);
            for i = 1:E.m
                Y(:,i) = E.Xi(1:6,i);
            end
            
            [y, Py] = unscented_transformation(E, Y, E.Wi, R);
            
            Pxz = zeros(E.n, 6);
            for i = 1:E.m
                Pxz = Pxz + E.Wi(i)*(E.Xi(:,i) - E.x)*(Y(:,i) - y)';
            end
            
            K = Pxz*inv(Py);
            
            E.x = E.x + K*(z - y);
            E.P = E.P - K*Py*K';
        end
        %% Equation of motion
        function Xk = f(E, X, h, a_imu, W_imu)
            xk = X(1:3);
            vk = X(4:6);
            ypr = X(7:9);
            Rk = eul2rotm(ypr');
            b_a_k = X(10);
            
            R_pre = Rk;
            b_a_pre = b_a_k;
            
            W = E.R_bi*W_imu;
            Rk = Rk*expm_SO3(h/2*(W + E.W_pre));
            
            % This assumes IMU provide acceleration without g
            ak = Rk*E.R_bi*a_imu + b_a_k*E.e3;
            a_pre = R_pre*E.R_bi*E.a_imu_pre + b_a_pre*E.e3;
               
            xk = xk + h*vk + h^2/2*a_pre;
            vk = vk + h/2*(ak + a_pre);
            
            E.W_pre = W;
            E.a_imu_pre = a_imu;
            E.a = ak;
          
            Xk = [xk; vk; rotm2eul(Rk)'; b_a_k];         
        end
        %% Delayed Prediction
        function delayed_prediction(E, h, a_imu, W_imu)
            
            E.delayed_sigma_points();
            
            Xk = zeros(E.n_s, E.m_s);
            for i = 1:E.m_s
                Xk(:,i) = delayed_f(E, E.Xi_s(:, i), h, a_imu, W_imu);
            end
            E.Xi_s = Xk;
            
            [E.x_s, E.P_s] = unscented_transformation(E, Xk, E.Wi_s, E.Q_s);
            E.x = E.x_s(1:E.n,1);
            E.P = E.P_s(1:E.n,1:E.n);
        end
        %% Delayed IMU correction
        function delayed_correction_imu(E, ypr, V_R_imu)
            
            % z = reshape(R, 9, 1);
            z = ypr;
            R = V_R_imu;
            
            Y = zeros(3, E.m_s);
            for i = 1:E.m_s
                Y(:,i) = E.Xi_s(7:9,i);
            end
            
            [y, Py] = unscented_transformation(E, Y, E.Wi_s, R);
            
            Pxz = zeros(E.n_s, 3);
            for i = 1:E.m_s
                Pxz = Pxz + E.Wi_s(i)*(E.Xi_s(:,i) - E.x_s)*(Y(:,i) - y)';
            end
            
            K = Pxz*inv(Py);
            
            E.x_s = E.x_s + K*(z - y);
            E.P_s = E.P_s - K*Py*K';
            
            E.x = E.x_s(1:E.n,1);
            E.P = E.P_s(1:E.n,1:E.n);
        end
        %% Delayed GPS correction
        function delayed_correction_gps(E, x_gps, v_gps, V_x_gps, V_v_gps)
            
            z = [x_gps; v_gps];
            R = blkdiag(V_x_gps, V_v_gps);
            
            Y = zeros(6, E.m_s);
            for i = 1:E.m_s
                Y(:,i) = E.Xi_s(11:16,i);
            end
            
            [y, Py] = unscented_transformation(E, Y, E.Wi_s, R);
            Pxz = E.P_s(1:10,11:16);

            K = Pxz*inv(Py);
            
            E.x = E.x_s(1:10) + K*(z - y);
            E.P = E.P_s(1:10,1:10) - K*Py*K';
        end
        %% Delayed equation of motion
        function Xk = delayed_f(E, X, h, a_imu, W_imu)
            xk = X(1:3);
            vk = X(4:6);
            ypr = real(X(7:9));
            Rk = eul2rotm(ypr');
            b_a_k = X(10);
            
            R_pre = Rk;
            b_a_pre = b_a_k;
            
            W = E.R_bi*W_imu;
            Rk = Rk*expm_SO3(h/2*(W + E.W_pre));
            
            % This assumes IMU provide acceleration without g
            ak = Rk*E.R_bi*a_imu + b_a_k*E.e3;
            a_pre = R_pre*E.R_bi*E.a_imu_pre + b_a_pre*E.e3;
               
            xk = xk + h*vk + h^2/2*a_pre;
            vk = vk + h/2*(ak + a_pre);
            
            % Save data
            E.W_pre = W;
            E.a_imu_pre = a_imu;
            E.a = ak;
          
            Xk = [xk; vk; rotm2eul(Rk)'; b_a_k; X(11:20)];         
        end
        %% Delayed sigma points
        function delayed_sigma_points(E)
            E.Xi_s(:,1) = E.x_s;
            E.Wi_s(1) = E.kappa / (E.n_s + E.kappa);
            
            sqrt_P = sqrtm(E.P_s);
            
            for i = 1:E.n_s
                E.Xi_s(:,i+1) = E.x_s + E.gamma*sqrt_P(:,i);
                E.Xi_s(:,i+E.n_s+1) = E.x_s - E.gamma*sqrt_P(:,i);
            end
            
            for i = 2:E.m_s
                E.Wi_s(i) = 1 / (2*(E.n_s+E.kappa));
                E.Wi_s(i+E.n_s+1) = 1 / (2*(E.n_s + E.kappa));
            end
        end
        %% Augment the state
        function update_augmented_state(E)
            E.x_s = [E.x; E.x];
            E.P_s = [E.P, E.P;
                E.P, E.P];
        end
        %% Sigma points
        function sigma_points(E)
            E.Xi(:,1) = E.x;
            E.Wi(1) = E.kappa / (E.n + E.kappa);
            
            sqrt_P = sqrtm(E.P);
            
            for i = 1:E.n
                E.Xi(:,i+1) = E.x + E.gamma*sqrt_P(:,i);
                E.Xi(:,i+E.n+1) = E.x - E.gamma*sqrt_P(:,i);
            end
            
            for i = 2:E.m
                E.Wi(i) = 1 / (2*(E.n+E.kappa));
                E.Wi(i+E.n+1) = 1 / (2*(E.n + E.kappa));
            end
        end
        %% Unscented transformation
        function [xUT, PUT] = unscented_transformation(E, XUT, WUT, Noise)
            
            [mUT, nUT] = size(XUT);
            
            xUT = zeros(mUT, 1);
            for i = 1:nUT
                xUT = xUT + WUT(i)*XUT(:,i);
            end
            
            PUT = zeros(mUT, mUT);
            for i = 1:nUT
                PUT = PUT + WUT(i)*(XUT(:,i) - xUT)*(XUT(:,i) - xUT)';
            end
            
            PUT = PUT + Noise;
        end
        %% Update parameters
        function update_parameters(E, R_bi)
            E.R_bi = R_bi;
        end
        %% Update initial values
        function update_init_values(E, P, Q)
            E.P = P;
            E.Q = Q;
            E.Q_s = blkdiag(Q, Q);
            E.P_s = [P, P;
                P, P];
        end
        %% Get states
        function [x, v, R, P] = output_states(E)
            x = E.x(1:3);
            v = E.x(4:6);
            R = eul2rotm(real(E.x(7:9))');
            P = E.P;
        end
    end
end
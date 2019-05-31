% RunSimulation.m
% This code requires the Multi-Parametric Toolbox
% (http://control.ee.ethz.ch/~mpt/).

clear;
load('data');

NSim = 100; % Number of time steps
NDir = 100; % Number of Directions
TD = zeros(1,NDir); % Target directions

% Neural network model
NW = 2500;                  % Number of neurons
W = zeros(NW, 10);          % Synaptic wights of hidden layer
W2 = (2*rand(6,NW)-1)./NW;  % Synaptic wights of output layer
r = zeros(NW, NDir);        % Neuronal activities

% Data strages
Kal = zeros(8, 4, NSim, NDir);    % Kalman gains
OutG = zeros(6, 4, NDir);           % Sensory feedback gains
X_r = zeros(10, NSim, NDir);        % Estimated states
X_r_ = X_r;                         % Actual States

for dir = 1:(NDir+1)
    TD(dir) = (dir-1)/NDir*2*pi; % Target torque direction
    xref = zeros(8,NSim);
    
    for k = 1:NSim
        target_tau = 1./(1 + exp(-1*(k-10)))*[cos(TD(dir)) sin(TD(dir))];
        xref(1:2,k) =target_tau';
    end       
   
    x = zeros(8, NSim);
    X = zeros(10, NSim);
    xhat = x;
    x_ = x;
    X_ = X;
    
    u = zeros(6,NSim);
    u_ = u;
    y = zeros(4,NSim);
    
    noise.u = 0.2;
    noise.t = 0.01;
    noise.f = 0.1;
          
    Yn = diag([noise.t noise.t noise.f noise.f]);
          
    K = zeros(8,4,NSim-1);
    L = zeros(6,10,NSim-1);
    L_ = L;
    P = zeros(8,8,NSim);
    P(1:4,1:4,1) = Yn*Yn';
    
    %% Running simulation ----------------------------------------
    for k = 1:NSim
        X(:,k) = [xhat(:,k); xref(1:2,k)];
        for i = 1:ctrl.optimizer.Num
            if AA{i}*X(:,k) <= bb{i}
                break;
            end
        end
        L(:,:,k) = -F{i}(1:6,:); % Feedback gain
        u(:,k) = -L(:,:,k)*X(:,k); % Motor command
        
        % Actial state
        X_(:,k) = [x_(:,k); xref(1:2,k)];
        for i = 1:ctrl.optimizer.Num
            if AA{i}*X_(:,k) <= bb{i}
                break;
            end
        end      
        L_(:,:,k) = -F{i}(1:6,:); % Feedback gain
        u_(:,k) = -L_(:,:,k)*X_(:,k); % Motor command
        
        for j = 1:6
            if u(j,k) < 0
                u(j,k) = 0;
            end
            if u_(j,k) < 0
                u_(j,k) = 0;
            end            
        end
        
        ctr_n = noise.u*diag(randn(1,6));
        U = (eye(6) + ctr_n)*u(:,k);      
        
        for j = 1:6
            if U(j,1) < 0
                U(j,1) = 0;
            end  
        end       
        
        x(:,k+1) = A*x(:,k) + B*U;
        y(:,k) = C*x(:,k) + Yn*randn(4,1);
        
        x_(:,k+1) = A*x_(:,k) + B*u_(:,k);
                
        % Kalman filter
        P_ = A*P(:,:,k)*A' + (B*noise.u*u(:,k))*(B*noise.u*u(:,k))';
        K(:,:,k) = P_*C'*pinv(C*P_*C' + Yn*Yn');
        P(:,:,k+1) = (eye(8) - K(:,:,k)*C)*P_;
        
        xhat(:,k+1) = A*xhat(:,k) + B*u(:,k) ...
            + K(:,:,k)*(y(:,k) - C*xhat(:,k)); 
        
    end
    %%  ----------------------------------------------------------
  
    % Store data
    Kal(:,:,:,dir) = K; % Kalman gain
    OutG(:,:,dir) = L(:,1:8,end)*K(:,:,end);% Output gain
    X_r(:,:,dir) = X;
    X_r_(:,:,dir) = X_;
        
    if dir <= NDir
        W = -pinv(W2)*L_(:,:,end);
        r(:,dir) = W*X_(:,end);    % neural activity
    end 
    
end

[~, nPD] = max(abs(r'));    % Neuronal PDs
[~, mPD_] = max(squeeze(X_r_(3:8,end,:)),[], 2);    % Muscle PDs
mPD = (mPD_-1)/NDir*2*pi;

%% plot ---------------------------------------------------------
time = 0:dt*1e+3:(NSim-1)*dt*1e+3;

% Torque profiles
dir = 13;
figure(1);

subplot(1,2,1);
plot(time, X_r(9,:,dir),'k:','LineWidth',0.1); hold on;
plot(time, X_r(1,:,dir),'-b');
plot(time, X_r(1,:,dir),'--r');
hold off; box off;
xlim([-100 600]); ylim([0 1.1]);
xlabel('time [ms]');ylabel('Torque [Nm]');
title('Shoulder-torque profile');

subplot(1,2,2);
plot(time,X_r(10,:,dir),'k:','LineWidth',0.1); hold on;
plot(time, X_r(2,:,dir),'-b');
plot(time, X_r(2,:,dir),'--r');
hold off; box off;
xlim([-100 600]); ylim([0 1.1]);
xlabel('time [ms]');
title('Elbow-torque profile');

% Kalman gains
figure(2);
subplot(2,2,1);
plot(time, squeeze(Kal(1,1,:,dir)),'-b'); hold on;
plot(time, squeeze(Kal(1,2,:,dir)),'--r'); 
hold off; box off;
xlim([-100 600]);
title('Kalman gains K11, K12');

subplot(2,2,3);
plot(time(1:NSim), squeeze(Kal(1,3,:,dir)),'-b'); hold on;
plot(time(1:NSim), squeeze(Kal(1,4,:,dir)),'--r'); 
hold off; box off;
xlim([-100 600]);
xlabel('time [ms]');
title('Kalman gains K13, K14');

subplot(2,2,2);
plot(time(1:NSim), squeeze(Kal(2,1,:,dir)),'-b'); hold on;
plot(time(1:NSim), squeeze(Kal(2,2,:,dir)),'--r'); 
hold off; box off;
xlim([-100 600]);
title('Kalman gains K21, K22');

subplot(2,2,4);
plot(time(1:NSim), squeeze(Kal(2,3,:,dir)),'-b'); hold on;
plot(time(1:NSim), squeeze(Kal(2,4,:,dir)),'--r'); 
hold off; box off;
xlim([-100 600]);
xlabel('time [ms]');
title('Kalman gains K23, K24');


% Muscle preferences
figure(3);
circle_x = 25*cos(TD); circle_y = 25*sin(TD);
range = 50;

subplot(2,3,1);
plot([0 0],[-range range],'k-', 'LineWidth',0.1); hold on;
plot([-range range],[0 0],'k-', 'LineWidth',0.1);
plot(circle_x,circle_y,'k-', 'LineWidth',0.1);
plot(squeeze(X_r_(3,end,:)).*(cos(TD))',  squeeze(X_r_(3,end,:)).*(sin(TD))','r.');
plot(squeeze(X_r_(3,end,:)).*(cos(TD))',  squeeze(X_r_(3,end,:)).*(sin(TD))','r','LineWidth',0.25);
plot(X_r_(3,end,mPD_(1)).*cos(mPD(1)),  X_r_(3,end,mPD_(1)).*sin(mPD(1)),'rx');
hold off; axis square; xlim([-range range]); ylim([-range range]);
box off; title('SF');

subplot(2,3,4);
plot([0 0],[-range range],'k-', 'LineWidth',0.1); hold on;
plot(1.2*[-range range],[0 0],'k-', 'LineWidth',0.1);
plot(circle_x,circle_y,'k-', 'LineWidth',0.1);
plot(squeeze(X_r_(4,end,:)).*(cos(TD))',  squeeze(X_r_(4,end,:)).*(sin(TD))','r.');
plot(squeeze(X_r_(4,end,:)).*(cos(TD))',  squeeze(X_r_(4,end,:)).*(sin(TD))','r','LineWidth',0.25);
plot(X_r_(4,end,mPD_(2)).*cos(mPD(2)),  X_r_(4,end,mPD_(2)).*sin(mPD(2)),'rx');
hold off; axis square;  xlim([-range range]); ylim([-range range]);
box off; title('SX');
xlabel('Shoulder torque');
ylabel('Elbow torque');

subplot(2,3,3);
plot([0 0],[-range range],'k-', 'LineWidth',0.1); hold on;
plot([-range range],[0 0],'k-', 'LineWidth',0.1);
plot(circle_x,circle_y,'k-', 'LineWidth',0.1);
plot(squeeze(X_r_(5,end,:)).*(cos(TD))',  squeeze(X_r_(5,end,:)).*(sin(TD))','b.');
plot(squeeze(X_r_(5,end,:)).*(cos(TD))',  squeeze(X_r_(5,end,:)).*(sin(TD))','b','LineWidth',0.25);
plot(X_r_(5, end,mPD_(3)).*cos(mPD(3)),  X_r_(5, end,mPD_(3)).*sin(mPD(3)),'bx');
hold off; axis square;  xlim([-range range]); ylim([-range range]);
box off; title('EF');

subplot(2,3,6);
plot([0 0],[-range range],'k-', 'LineWidth',0.1); hold on;
plot([-range range],[0 0],'k-', 'LineWidth',0.1);
plot(circle_x,circle_y,'k-', 'LineWidth',0.1);
plot(squeeze(X_r_(6,end,:)).*(cos(TD))',  squeeze(X_r_(6,end,:)).*(sin(TD))','b.');
plot(squeeze(X_r_(6,end,:)).*(cos(TD))',  squeeze(X_r_(6,end,:)).*(sin(TD))','b','LineWidth',0.25);
plot(X_r_(6,end,mPD_(4)).*cos(mPD(4)),  X_r_(6,end,mPD_(4)).*sin(mPD(4)),'bx');
hold off; axis square; xlim([-range range]); ylim([-range range]);
box off; title('EX');

subplot(2,3,2);
plot([0 0],[-range range],'k-', 'LineWidth',0.1); hold on;
plot([-range range],[0 0],'k-', 'LineWidth',0.1);
plot(circle_x,circle_y,'k-', 'LineWidth',0.1);
plot(squeeze(X_r_(7,end,:)).*(cos(TD))',  squeeze(X_r_(7,end,:)).*(sin(TD))','g.');
plot(squeeze(X_r_(7,end,:)).*(cos(TD))',  squeeze(X_r_(7,end,:)).*(sin(TD))','g','LineWidth',0.25);
plot(X_r_(7,end,mPD_(5)).*cos(mPD(5)),  X_r_(7,end,mPD_(5)).*sin(mPD(5)),'gx');
hold off; axis square;  xlim([-range range]); ylim([-range range]);
box off; title('BF');

subplot(2,3,5);
plot([0 0],[-range range],'k-', 'LineWidth',0.1); hold on;
plot([-range range],[0 0],'k-', 'LineWidth',0.1);
plot(circle_x,circle_y,'k-', 'LineWidth',0.1);
plot(squeeze(X_r_(8,end,:)).*(cos(TD))',  squeeze(X_r_(8,end,:)).*(sin(TD))','g.');
plot(squeeze(X_r_(8,end,:)).*(cos(TD))',  squeeze(X_r_(8,end,:)).*(sin(TD))','g','LineWidth',0.25);
plot(X_r_(8,end,mPD_(6)).*cos(mPD(6)),  X_r_(8,end,mPD_(6)).*sin(mPD(6)),'gx');
hold off; axis square;  xlim([-range range]); ylim([-range range]);
box off; title('BX');


% Distribution of neuronal PDs
bin = (0:19)/20*2*pi;
figure(4);
subplot(1,2,1);
plot(W(:,1),W(:,2),'r.', 'MarkerSize',4); axis image;box off;
xlabel('Shoulder torque');ylabel('Elbow torque');
title('Synaptic weights');

subplot(1,2,2);
rose((nPD-1)*2*pi/NDir, bin);
title('Distribution of neuronal PDs');

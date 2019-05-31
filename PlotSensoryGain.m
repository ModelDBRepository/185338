% PlotSensoryGain.m
% This code requires the Multi-Parametric Toolbox
% (http://control.ee.ethz.ch/~mpt/).

clear;
load('data');

NDir = 16;	% Number of Directions
NSim = 100; % Number of time steps

HandDir = zeros(1, NDir+1); % Desired hand direction
Fd = zeros(2, NDir+1);      % Desired hand force
Td = zeros(2, NDir+1);      % Desired torque
Gain = zeros(6, 4, NDir+1); % Sensory feedback gain
Gain_ = Gain;

PV = zeros(6, 2, NSim, NDir);
PV_ = zeros(6, 2, NSim, NDir);
Fdyn = zeros(2, NSim+1, NDir);
Fdyn_ = zeros(2, NSim+1, NDir);

for dir = 1:NDir
    % Desired motor output
    HandDir(dir) = (dir-1)/NDir*2*pi;
    xref = zeros(8,NSim);
    
    % Desired force
    Fd(:,dir) = 3*[cos(HandDir(dir)); sin(HandDir(dir))];
    % Desired torque
    Td(:,dir) = (Jacob'*Fd(:,dir));          
        
   %% initialize policy and filter
    x = zeros(8, NSim);
    xhat = x;
     
    noise.u = 0.2;
    noise.t = 0.01;
    noise.f = 0.1;
    noise.w = 0.2;
    
    Cn = noise.u;
    Yn = diag([noise.t noise.t noise.f noise.f]); 
    Wn = noise.w;
    
    K = zeros(8, 4, NSim-1);
    L = zeros(6, 10, NSim-1);
    
    P = zeros(8,8);
    P(1:4,1:4) = Yn*Yn';
    
    % For extra condition without SDN
    x_ = x; xhat_ = x_;
    K_ = K; L_ = L; P_ = P;
    
    %% Running Simulation
    for k = 1:NSim
        % Desired torque
        xref(1:2,k) = 1./(1 + exp(-0.5*(k-10)))*Td(:,dir);
                
       %% Normal condition with SDN
        X = [xhat(:,k); xref(1:2,k)];
        
        for i = 1:ctrl.optimizer.Num
            if AA{i}*X <= bb{i}
                break;
            end
        end
        
        L(:,:,k) = -F{i}(1:6,:);        
        u = -L(:,:,k)*X;
        
        for j = 1:6
            if u(j,1) < 0
                u(j,1) = 0;
            end
        end
        
        U = (eye(6) + Cn*diag(randn(1,6)))*u;
        for j = 1:6
            if U(j,1) < 0
                U(j,1) = 0;
            end
        end
        
        x(:,k+1) = A*x(:,k) + B*U;
        y = C*x(:,k) + Yn*randn(4,1);
        
        % kalman filter
        P = A*P*A' + (B*Cn*u)*(B*Cn*u)';
        K(:,:,k) = P*C'*pinv(C*P*C' + Yn*Yn');
        P = (eye(8) - K(:,:,k)*C)*P;
        
        xhat(:,k+1) = A*xhat(:,k) + B*u + K(:,:,k)*(y - C*xhat(:,k));
        
       %% Extra condition with non-SDN
        X = [xhat_(:,k); xref(1:2,k)];
        
        for i = 1:ctrl.optimizer.Num
            if AA{i}*X <= bb{i}
                break;
            end
        end
        
        L_(:,:,k) = -F{i}(1:6,:);        
        u = -L_(:,:,k)*X;
        
        for j = 1:6
            if u(j,1) < 0
                u(j,1) = 0;
            end
        end
        
        U = u + Wn*randn(6,1);
        for j = 1:6
            if U(j,1) < 0
                U(j,1) = 0;
            end
        end
        
        x_(:,k+1) = A*x_(:,k) + B*U;
        y_ = C*x_(:,k) + Yn*randn(4,1);
        
        % kalman filter
        P_ = A*P_*A' + (B*Wn)*(B*Wn)';
        K_(:,:,k) = P_*C'*pinv(C*P_*C' + Yn*Yn');
        P_ = (eye(8) - K_(:,:,k)*C)*P_;
        
        xhat_(:,k+1) = A*xhat_(:,k) + B*u + K_(:,:,k)*(y_ - C*xhat_(:,k));
        
        % Population vector
        PV(:,:,k,dir) = L(:,1:8,k)*K(:,3:4,k);
        PV_(:,:,k,dir) = L_(:,1:8,k)*K_(:,3:4,k);       
    end
    
    % output gain    
    Gain(:,:,dir) = L(:,1:8,end)*K(:,:,end);
    Gain_(:,:,dir) = L_(:,1:8,end)*K_(:,:,end);
    
    % Force dynamics
    Fdyn(:,:,dir) = pinv(Jacob')*x(1:2,:);
    Fdyn_(:,:,dir) = pinv(Jacob')*x_(1:2,:);    
    
end

HandDir(:,end) = 2*pi;
Fd(:,end) = Fd(:,1);
Td(:,end) = Td(:,1);
Gain(:,:,end) = Gain(:,:,1);

%% Plot figures of normal condtion with SDN ----------------------
% Distributions of the sensory feedback gains
figure(1);
subplot(1,2,1);
tmp = [1 4 10];
hold on;
for k = tmp
    plot(1000*[0 Td(1,k)],1000*[0 Td(2, k)],':k');
end
plot(squeeze(Gain(1,1,tmp)),squeeze(Gain(1,2,tmp)),'ro');
plot(squeeze(Gain(2,1,tmp)),squeeze(Gain(2,2,tmp)),'rx');
plot(squeeze(Gain(3,1,tmp)),squeeze(Gain(3,2,tmp)),'go');
plot(squeeze(Gain(4,1,tmp)),squeeze(Gain(4,2,tmp)),'gx');
plot(squeeze(Gain(5,1,tmp)),squeeze(Gain(5,2,tmp)),'bo');
plot(squeeze(Gain(6,1,tmp)),squeeze(Gain(6,2,tmp)),'bx');

title('Torque gain with SDN');
xlim([-1000 1000]);
ylim([-1000 1000]);
xlabel('Shoulder torque'); 
ylabel('Elbow torque');
axis square;
hold off;

subplot(1,2,2);
hold on;
for k = tmp
    plot(5*[0 Fd(1,k)],5*[0 Fd(2,k)],':k');
end
plot(squeeze(Gain(1,3,tmp)),squeeze(Gain(1,4,tmp)),'ro');
plot(squeeze(Gain(2,3,tmp)),squeeze(Gain(2,4,tmp)),'rx');
plot(squeeze(Gain(3,3,tmp)),squeeze(Gain(3,4,tmp)),'go');
plot(squeeze(Gain(4,3,tmp)),squeeze(Gain(4,4,tmp)),'gx');
plot(squeeze(Gain(5,3,tmp)),squeeze(Gain(5,4,tmp)),'bo');
plot(squeeze(Gain(6,3,tmp)),squeeze(Gain(6,4,tmp)),'bx');

title('Force gain with SDN');
xlim([-20 20]);
ylim([-20 20]);
xlabel('x-force');
ylabel('y-force');
axis square;
hold off;

% Amplitudes of the sensory feedback gains
figure(2);
subplot(2,3,1);
% Below commenting out code requires MATLAB Curve Fitting Toolbox.
% clear x y;
% j = 1;
% for k=1:dir
%     tmp = squeeze(sqrt(Gain(1,3,k).^2+Gain(1,4,k).^2));
%     if tmp > 0.001
%         x(j,1) = HandDir(k);
%         y(j,1) = tmp;
%         j = j+1;
%     end
% end
% type = fittype('a*cos(x+c)');
% mycurve = fit(x,y,type);
% plot([0 360], zeros(2), 'k:',...
%     HandDir*180/pi, mycurve(HandDir'), 'm--',...
%     HandDir*180/pi, squeeze(sqrt(Gain(1,3,:).^2+Gain(1,4,:).^2)), 'rx');
plot([0 360], zeros(2), 'k:',...
     HandDir*180/pi, squeeze(sqrt(Gain(1,3,:).^2+Gain(1,4,:).^2)), 'rx');
box off;
xlim([0 360]); ylim([-20 20]);
set(gca,'xtick',0:180:360);
title('Gain amplitudes with SDN');

subplot(2,3,2);
% Below commenting out code requires MATLAB Curve Fitting Toolbox.
% clear x y;
% j = 1;
% for k=1:dir
%     tmp = squeeze(sqrt(Gain(3,3,k).^2+Gain(3,4,k).^2));
%     if tmp > 0.001
%         x(j,1) = HandDir(k);
%         y(j,1) = tmp;
%         j = j+1;
%     end
% end
% mycurve = fit(x,y,type);
% plot([0 360], zeros(2), 'k:',...
%     HandDir*180/pi, mycurve(HandDir'), 'm--',...
%     HandDir*180/pi, squeeze(sqrt(Gain(3,3,:).^2+Gain(3,4,:).^2)), 'gx');
 plot([0 360], zeros(2), 'k:',...
     HandDir*180/pi, squeeze(sqrt(Gain(3,3,:).^2+Gain(3,4,:).^2)), 'gx');
box off;
xlim([0 360]); ylim([-20 20]);
set(gca,'xtick',0:180:360);

subplot(2,3,3);
% Below commenting out code requires MATLAB Curve Fitting Toolbox.
% clear x y;
% j = 1;
% for k=1:dir
%     tmp = squeeze(sqrt(Gain(5,3,k).^2+Gain(5,4,k).^2));
%     if tmp > 0.001
%         x(j,1) = HandDir(k);
%         y(j,1) = tmp;
%         j = j+1;
%     end
% end
% mycurve = fit(x,y,type);
% plot([0 360], zeros(2), 'k:',...
%     HandDir*180/pi, mycurve(HandDir'), 'm--',...
%     HandDir*180/pi, squeeze(sqrt(Gain(5,3,:).^2+Gain(5,4,:).^2)), 'bx');
plot([0 360], zeros(2), 'k:',...
     HandDir*180/pi, squeeze(sqrt(Gain(5,3,:).^2+Gain(5,4,:).^2)), 'bx');
box off;
xlim([0 360]); ylim([-20 20]);
set(gca,'xtick',0:180:360);

subplot(2,3,4);
% Below commenting out code requires MATLAB Curve Fitting Toolbox.
% clear x y;
% j = 1;
% for k=1:dir
%     tmp = squeeze(sqrt(Gain(2,3,k).^2+Gain(2,4,k).^2));
%     if tmp > 0.001
%         x(j,1) = HandDir(k);
%         y(j,1) = tmp;
%         j = j+1;
%     end
% end
% mycurve = fit(x,y,type);
% plot([0 360], zeros(2), 'k:',...
%     HandDir*180/pi, mycurve(HandDir'), 'm--',...
%     HandDir*180/pi, squeeze(sqrt(Gain(2,3,:).^2+Gain(2,4,:).^2)), 'rx');
plot([0 360], zeros(2), 'k:',...
     HandDir*180/pi, squeeze(sqrt(Gain(2,3,:).^2+Gain(2,4,:).^2)), 'rx');
box off;
xlim([0 360]); ylim([-20 20]);
set(gca,'xtick',0:180:360);
xlabel('Force direction [degs]');
ylabel('Amplitude of the gain');

subplot(2,3,5);
% Below commenting out code requires MATLAB Curve Fitting Toolbox.
% clear x y;
% j = 1;
% for k=1:dir
%     tmp = squeeze(sqrt(Gain(4,3,k).^2+Gain(4,4,k).^2));
%     if tmp > 0.001
%         x(j,1) = HandDir(k);
%         y(j,1) = tmp;
%         j = j+1;
%     end
% end
% mycurve = fit(x,y,type);
% plot([0 360], zeros(2), 'k:',...
%     HandDir*180/pi, mycurve(HandDir'), 'm--',...
%     HandDir*180/pi, squeeze(sqrt(Gain(4,3,:).^2+Gain(4,4,:).^2)), 'gx');
plot([0 360], zeros(2), 'k:',...
     HandDir*180/pi, squeeze(sqrt(Gain(4,3,:).^2+Gain(4,4,:).^2)), 'gx');
box off;
xlim([0 360]); ylim([-20 20]);
set(gca,'xtick',0:180:360);

subplot(2,3,6);
% Below commenting out code requires MATLAB Curve Fitting Toolbox.
% clear x y;
% j = 1;
% for k=1:dir
%     tmp = squeeze(sqrt(Gain(6,3,k).^2+Gain(6,4,k).^2));
%     if tmp > 0.001
%         x(j,1) = HandDir(k);
%         y(j,1) = tmp;
%         j = j+1;
%     end
% end
% mycurve = fit(x,y,type);
% plot([0 360], zeros(2), 'k:',...
%     HandDir*180/pi, mycurve(HandDir'), 'm--',...
%     HandDir*180/pi, squeeze(sqrt(Gain(6,3,:).^2+Gain(6,4,:).^2)), 'bx');
plot([0 360], zeros(2), 'k:',...
     HandDir*180/pi, squeeze(sqrt(Gain(6,3,:).^2+Gain(6,4,:).^2)), 'bx');
box off;
xlim([0 360]); ylim([-20 20]);
set(gca,'xtick',0:180:360);

% Force gain preferences
figure(3);
amp = max(max(sqrt(Gain(1:6,3,:).^2+Gain(1:6,4,:).^2)));
subplot(2,3,1);
plot([0 amp*cos(0)], [0, amp*sin(0)],'k:'); hold on;
for k=1:3
    plot([0 amp*cos((k)/2*pi)], [0, amp*sin((k)/2*pi)],'k:');
end
plot(amp*cos((0:100)./100*2*pi), amp*sin((0:100)./100*2*pi),'k-');
plot(squeeze(Gain(1,3,:))',squeeze(Gain(1,4,:))','or-');
hold off; box off; axis square;
xlim([-20 20]); ylim([-20 20]);
title('Gain preferences with SDN');

subplot(2,3,4);
plot([0 amp*cos(0)], [0, amp*sin(0)],'k:'); hold on;
for k=1:3
    plot([0 amp*cos((k)/2*pi)], [0, amp*sin((k)/2*pi)],'k:');
end
plot(amp*cos((0:100)./100*2*pi), amp*sin((0:100)./100*2*pi),'k-');
plot(squeeze(Gain(2,3,:))',squeeze(Gain(2,4,:))','xr-');
hold off; box off; axis square;
xlim([-20 20]); ylim([-20 20]);
xlabel('x-force gain'); ylabel('y-force gain');

subplot(2,3,3);
plot([0 amp*cos(0)], [0, amp*sin(0)],'k:'); hold on;
for k=1:3
    plot([0 amp*cos((k)/2*pi)], [0, amp*sin((k)/2*pi)],'k:');
end
plot(amp*cos((0:100)./100*2*pi), amp*sin((0:100)./100*2*pi),'k-');
plot(squeeze(Gain(3,3,:))',squeeze(Gain(3,4,:))','ob-');
hold off; box off; axis square;
xlim([-20 20]); ylim([-20 20]);

subplot(2,3,6);
plot([0 amp*cos(0)], [0, amp*sin(0)],'k:'); hold on;
for k=1:3
    plot([0 amp*cos((k)/2*pi)], [0, amp*sin((k)/2*pi)],'k:');
end
plot(amp*cos((0:100)./100*2*pi), amp*sin((0:100)./100*2*pi),'k-');
plot(squeeze(Gain(4,3,:))',squeeze(Gain(4,4,:))','xb-');
hold off; box off; axis square;
xlim([-20 20]); ylim([-20 20]);

subplot(2,3,2);
plot([0 amp*cos(0)], [0, amp*sin(0)],'k:'); hold on;
for k=1:3
    plot([0 amp*cos((k)/2*pi)], [0, amp*sin((k)/2*pi)],'k:');
end
plot(amp*cos((0:100)./100*2*pi), amp*sin((0:100)./100*2*pi),'k-');
plot(squeeze(Gain(5,3,:))',squeeze(Gain(5,4,:))','og-');
hold off; box off; axis square;
xlim([-20 20]); ylim([-20 20]);

subplot(2,3,5);
plot([0 amp*cos(0)], [0, amp*sin(0)],'k:'); hold on;
for k=1:3
    plot([0 amp*cos((k)/2*pi)], [0, amp*sin((k)/2*pi)],'k:');
end
plot(amp*cos((0:100)./100*2*pi), amp*sin((0:100)./100*2*pi),'k-');
plot(squeeze(Gain(6,3,:))',squeeze(Gain(6,4,:))','xg-');
hold off; box off; axis square;
xlim([-20 20]); ylim([-20 20]);


%% Plot figures of extra condtion with non-SDN -------------------
% Distributions of the sensory feedback gains
figure(11);
subplot(1,2,1);
tmp = [1 4 10];
hold on;
for k = tmp
    plot(1000*[0 Td(1,k)],1000*[0 Td(2, k)],':k');
end
plot(squeeze(Gain_(1,1,tmp)),squeeze(Gain_(1,2,tmp)),'ro');
plot(squeeze(Gain_(2,1,tmp)),squeeze(Gain_(2,2,tmp)),'rx');
plot(squeeze(Gain_(3,1,tmp)),squeeze(Gain_(3,2,tmp)),'go');
plot(squeeze(Gain_(4,1,tmp)),squeeze(Gain_(4,2,tmp)),'gx');
plot(squeeze(Gain_(5,1,tmp)),squeeze(Gain_(5,2,tmp)),'bo');
plot(squeeze(Gain_(6,1,tmp)),squeeze(Gain_(6,2,tmp)),'bx');

title('Torque gain with non-SDN');
xlim([-50 50]);
ylim([-50 50]);
xlabel('Shoulder torque'); 
ylabel('Elbow torque');
axis square;
hold off;

subplot(1,2,2);
hold on;
for k = tmp
    plot(5*[0 Fd(1,k)],5*[0 Fd(2,k)],':k');
end
plot(squeeze(Gain_(1,3,tmp)),squeeze(Gain_(1,4,tmp)),'ro');
plot(squeeze(Gain_(2,3,tmp)),squeeze(Gain_(2,4,tmp)),'rx');
plot(squeeze(Gain_(3,3,tmp)),squeeze(Gain_(3,4,tmp)),'go');
plot(squeeze(Gain_(4,3,tmp)),squeeze(Gain_(4,4,tmp)),'gx');
plot(squeeze(Gain_(5,3,tmp)),squeeze(Gain_(5,4,tmp)),'bo');
plot(squeeze(Gain_(6,3,tmp)),squeeze(Gain_(6,4,tmp)),'bx');

title('Force gain with non-SDN');
xlim([-1 1]);
ylim([-1 1]);
xlabel('x-force');
ylabel('y-force');
axis square;
hold off;

% Amplitudes of the sensory feedback gains
figure(12);
subplot(2,3,1);

% Below commenting out code requires MATLAB Curve Fitting Toolbox.
% clear x y;
% j = 1;
% for k=1:dir
%     tmp = squeeze(sqrt(Gain_(1,3,k).^2+Gain_(1,4,k).^2));
%     if tmp > 0.001
%         x(j,1) = HandDir(k);
%         y(j,1) = tmp;
%         j = j+1;
%     end
% end
% mycurve = fit(x,y,type);
% plot([0 360], zeros(2), 'k:',...
%     HandDir*180/pi, mycurve(HandDir'), 'm--',...
%     HandDir*180/pi, squeeze(sqrt(Gain_(1,3,:).^2+Gain_(1,4,:).^2)), 'rx');
plot([0 360], zeros(2), 'k:',...    
    HandDir*180/pi, squeeze(sqrt(Gain_(1,3,:).^2+Gain_(1,4,:).^2)), 'rx');
box off;
xlim([0 360]); ylim([-1 1]);
set(gca,'xtick',0:180:360);
title('Gain amplitudes with non-SDN');

subplot(2,3,2);
% Below commenting out code requires MATLAB Curve Fitting Toolbox.
% clear x y;
% j = 1;
% for k=1:dir
%     tmp = squeeze(sqrt(Gain_(3,3,k).^2+Gain_(3,4,k).^2));
%     if tmp > 0.001
%         x(j,1) = HandDir(k);
%         y(j,1) = tmp;
%         j = j+1;
%     end
% end
% mycurve = fit(x,y,type);
% plot([0 360], zeros(2), 'k:',...
%     HandDir*180/pi, mycurve(HandDir'), 'm--',...
%     HandDir*180/pi, squeeze(sqrt(Gain_(3,3,:).^2+Gain_(3,4,:).^2)), 'gx')
plot([0 360], zeros(2), 'k:',...
    HandDir*180/pi, squeeze(sqrt(Gain_(3,3,:).^2+Gain_(3,4,:).^2)), 'gx')
box off;
xlim([0 360]); ylim([-1 1]);
set(gca,'xtick',0:180:360);

subplot(2,3,3);
% Below commenting out code requires MATLAB Curve Fitting Toolbox.
% clear x y;
% j = 1;
% for k=1:dir
%     tmp = squeeze(sqrt(Gain_(5,3,k).^2+Gain_(5,4,k).^2));
%     if tmp > 0.001
%         x(j,1) = HandDir(k);
%         y(j,1) = tmp;
%         j = j+1;
%     end
% end
% mycurve = fit(x,y,type);
% plot([0 360], zeros(2), 'k:',...
%     HandDir*180/pi, mycurve(HandDir'), 'm--',...
%     HandDir*180/pi, squeeze(sqrt(Gain_(5,3,:).^2+Gain_(5,4,:).^2)), 'bx');
plot([0 360], zeros(2), 'k:',...
    HandDir*180/pi, squeeze(sqrt(Gain_(5,3,:).^2+Gain_(5,4,:).^2)), 'bx');
box off;
xlim([0 360]); ylim([-1 1]);
set(gca,'xtick',0:180:360);

subplot(2,3,4);
% Below commenting out code requires MATLAB Curve Fitting Toolbox.
% clear x y;
% j = 1;
% for k=1:dir
%     tmp = squeeze(sqrt(Gain_(2,3,k).^2+Gain_(2,4,k).^2));
%     if tmp > 0.001
%         x(j,1) = HandDir(k);
%         y(j,1) = tmp;
%         j = j+1;
%     end
% end
% mycurve = fit(x,y,type);
% plot([0 360], zeros(2), 'k:',...
%     HandDir*180/pi, mycurve(HandDir'), 'm--',...
%     HandDir*180/pi, squeeze(sqrt(Gain_(2,3,:).^2+Gain_(2,4,:).^2)), 'rx');
plot([0 360], zeros(2), 'k:',...
    HandDir*180/pi, squeeze(sqrt(Gain_(2,3,:).^2+Gain_(2,4,:).^2)), 'rx');
box off;
xlim([0 360]); ylim([-1 1]);
set(gca,'xtick',0:180:360);
xlabel('Force direction [degs]');
ylabel('Amplitude of the gain');

subplot(2,3,5);
% Below commenting out code requires MATLAB Curve Fitting Toolbox.
% clear x y;
% j = 1;
% for k=1:dir
%     tmp = squeeze(sqrt(Gain_(4,3,k).^2+Gain_(4,4,k).^2));
%     if tmp > 0.001
%         x(j,1) = HandDir(k);
%         y(j,1) = tmp;
%         j = j+1;
%     end
% end
% mycurve = fit(x,y,type);
% plot([0 360], zeros(2), 'k:',...
%     HandDir*180/pi, mycurve(HandDir'), 'm--',...
%     HandDir*180/pi, squeeze(sqrt(Gain_(4,3,:).^2+Gain_(4,4,:).^2)), 'gx');
plot([0 360], zeros(2), 'k:',...
     HandDir*180/pi, squeeze(sqrt(Gain_(4,3,:).^2+Gain_(4,4,:).^2)), 'gx');
box off;
box off;
xlim([0 360]); ylim([-1 1]);
set(gca,'xtick',0:180:360);

subplot(2,3,6);
% Below commenting out code requires MATLAB Curve Fitting Toolbox.
% clear x y;
% j = 1;
% for k=1:dir
%     tmp = squeeze(sqrt(Gain_(6,3,k).^2+Gain_(6,4,k).^2));
%     if tmp > 0.001
%         x(j,1) = HandDir(k);
%         y(j,1) = tmp;
%         j = j+1;
%     end
% end
% mycurve = fit(x,y,type);
% plot([0 360], zeros(2), 'k:',...
%     HandDir*180/pi, mycurve(HandDir'), 'm--',...
%     HandDir*180/pi, squeeze(sqrt(Gain_(6,3,:).^2+Gain_(6,4,:).^2)), 'bx');
plot([0 360], zeros(2), 'k:',...
    HandDir*180/pi, squeeze(sqrt(Gain_(6,3,:).^2+Gain_(6,4,:).^2)), 'bx');
box off;
xlim([0 360]); ylim([-1 1]);
set(gca,'xtick',0:180:360);

% Force Gain preferences
figure(13);
amp = max(max(sqrt(Gain_(1:6,3,:).^2+Gain_(1:6,4,:).^2)));
subplot(2,3,1);
plot([0 amp*cos(0)], [0, amp*sin(0)],'k:'); hold on;
for k=1:3
    plot([0 amp*cos((k)/2*pi)], [0, amp*sin((k)/2*pi)],'k:');
end
plot(amp*cos((0:100)./100*2*pi), amp*sin((0:100)./100*2*pi),'k-');
plot(squeeze(Gain_(1,3,:))',squeeze(Gain_(1,4,:))','or-');
hold off; box off; axis square;
xlim([-1 1]); ylim([-1 1]);
title('Gain preferences with non-SDN');

subplot(2,3,4);
plot([0 amp*cos(0)], [0, amp*sin(0)],'k:'); hold on;
for k=1:3
    plot([0 amp*cos((k)/2*pi)], [0, amp*sin((k)/2*pi)],'k:');
end
plot(amp*cos((0:100)./100*2*pi), amp*sin((0:100)./100*2*pi),'k-');
plot(squeeze(Gain_(2,3,:))',squeeze(Gain_(2,4,:))','xr-');
hold off; box off; axis square;
xlim([-1 1]); ylim([-1 1]);
xlabel('x-force gain'); ylabel('y-force gain');

subplot(2,3,3);
plot([0 amp*cos(0)], [0, amp*sin(0)],'k:'); hold on;
for k=1:3
    plot([0 amp*cos((k)/2*pi)], [0, amp*sin((k)/2*pi)],'k:');
end
plot(amp*cos((0:100)./100*2*pi), amp*sin((0:100)./100*2*pi),'k-');
plot(squeeze(Gain_(3,3,:))',squeeze(Gain_(3,4,:))','ob-');
hold off; box off; axis square;
xlim([-1 1]); ylim([-1 1]);

subplot(2,3,6);
plot([0 amp*cos(0)], [0, amp*sin(0)],'k:'); hold on;
for k=1:3
    plot([0 amp*cos((k)/2*pi)], [0, amp*sin((k)/2*pi)],'k:');
end
plot(amp*cos((0:100)./100*2*pi), amp*sin((0:100)./100*2*pi),'k-');
plot(squeeze(Gain_(4,3,:))',squeeze(Gain_(4,4,:))','xb-');
hold off; box off; axis square;
xlim([-1 1]); ylim([-1 1]);

subplot(2,3,2);
plot([0 amp*cos(0)], [0, amp*sin(0)],'k:'); hold on;
for k=1:3
    plot([0 amp*cos((k)/2*pi)], [0, amp*sin((k)/2*pi)],'k:');
end
plot(amp*cos((0:100)./100*2*pi), amp*sin((0:100)./100*2*pi),'k-');
plot(squeeze(Gain_(5,3,:))',squeeze(Gain_(5,4,:))','og-');
hold off; box off; axis square;
xlim([-1 1]); ylim([-1 1]);

subplot(2,3,5);
plot([0 amp*cos(0)], [0, amp*sin(0)],'k:'); hold on;
for k=1:3
    plot([0 amp*cos((k)/2*pi)], [0, amp*sin((k)/2*pi)],'k:');
end
plot(amp*cos((0:100)./100*2*pi), amp*sin((0:100)./100*2*pi),'k-');
plot(squeeze(Gain_(6,3,:))',squeeze(Gain_(6,4,:))','xg-');
hold off; box off; axis square;
xlim([-1 1]); ylim([-1 1]);


%% Muscle mechanical direction (Moment arms)
figure(4);
MD = Jacob'\J;
hand_pos = [l1_*cos(theta1)+l2_*cos(theta1+theta2), ...
    l1_*sin(theta1)+l2_*sin(theta1+theta2)];
subplot(1,4,1);
% Arm posture
plot([-hand_pos(1), l1_*cos(theta1)-hand_pos(1), 0], ...
    [-hand_pos(2),l1_*sin(theta1)-hand_pos(2), 0], 'k-'); 
hold on;
plot([0, MD(1,1)], [0,MD(2,1)], 'r-'); 
plot([0, MD(1,2)], [0,MD(2,2)], 'r--');
plot([0, MD(1,5)], [0,MD(2,5)], 'g-');
plot([0, MD(1,6)], [0,MD(2,6)], 'g--');
plot([0, MD(1,3)], [0,MD(2,3)], 'b-');
plot([0, MD(1,4)], [0,MD(2,4)], 'b--'); hold off;
axis square;
box off;
xlim([-0.25 0.25]); ylim([-0.25 0.25]);
xlabel('x-posiiton [m]');
ylabel('y-posiiton [m]');

subplot(1,4,2);
% Arm posture
plot([-hand_pos(1), l1_*cos(theta1)-hand_pos(1), 0], ...
    [-hand_pos(2),l1_*sin(theta1)-hand_pos(2), 0], 'k-'); 
hold on;
plot([0, 0.05*mean(Gain(1,3,:))], [0, 0.05*mean(Gain(1,4,:))], 'r-'); 
plot([0, 0.05*mean(Gain(2,3,:))], [0, 0.05*mean(Gain(2,4,:))], 'r--');
plot([0, 0.05*mean(Gain(5,3,:))], [0, 0.05*mean(Gain(5,4,:))], 'g-');
plot([0, 0.05*mean(Gain(6,3,:))], [0, 0.05*mean(Gain(6,4,:))], 'g--');
plot([0, 0.05*mean(Gain(3,3,:))], [0, 0.05*mean(Gain(3,4,:))], 'b-');
plot([0, 0.05*mean(Gain(4,3,:))], [0, 0.05*mean(Gain(4,4,:))], 'b--'); hold off;
axis square;
box off;
xlim([-0.25 0.25]); ylim([-0.25 0.25]);

subplot(1,4,3);
% Arm posture
plot([-hand_pos(1), l1_*cos(theta1)-hand_pos(1), 0], ...
    [-hand_pos(2),l1_*sin(theta1)-hand_pos(2), 0], 'k-'); 
hold on;
plot([0, 0.5*mean(Gain_(1,3,:))], [0, 0.5*mean(Gain_(1,4,:))], 'r-'); 
plot([0, 0.5*mean(Gain_(2,3,:))], [0, 0.5*mean(Gain_(2,4,:))], 'r--');
plot([0, 0.5*mean(Gain_(5,3,:))], [0, 0.5*mean(Gain_(5,4,:))], 'g-');
plot([0, 0.5*mean(Gain_(6,3,:))], [0, 0.5*mean(Gain_(6,4,:))], 'g--');
plot([0, 0.5*mean(Gain_(3,3,:))], [0, 0.5*mean(Gain_(3,4,:))], 'b-');
plot([0, 0.5*mean(Gain_(4,3,:))], [0, 0.5*mean(Gain_(4,4,:))], 'b--'); hold off;
axis square;
box off;
xlim([-0.25 0.25]); ylim([-0.25 0.25]);

subplot(1, 4, 4);
MD_deg =  atan2( MD(2,:), MD(1,:))'*180/pi;
Gain_deg =  atan2(mean(Gain(:,4,:),3), mean(Gain(:,3,:),3))*180/pi;
Gain_deg_ =  atan2(mean(Gain_(:,4,:),3), mean(Gain_(:,3,:),3))*180/pi;
errorbar(1:2, [mean(abs(Gain_deg-MD_deg)), mean(abs(Gain_deg_-MD_deg))],...
    [std(abs(Gain_deg-MD_deg)), std(abs(Gain_deg_-MD_deg))], 'b.-');
xlim([0.5 2.5]);
ylabel('Absolute angular error [deg]');
box off;
axis square;

% Force vector
figure(5);
Nv = 20;
subplot(3,2,1);
quiver(0:10:10*Nv, zeros(1,Nv+1),squeeze(Fdyn(1,1:Nv+1,7)),squeeze(Fdyn(2,1:Nv+1,7)));
box off;
xlim([0 200]);
ylim([-10 30]);
subplot(3,2,2);
quiver(0:10:10*Nv, zeros(1,Nv+1),squeeze(Fdyn_(1,1:Nv+1,7)),squeeze(Fdyn_(2,1:Nv+1,7)));
box off;
xlim([0 200]);
ylim([-10 30]);
subplot(3,2,3);
quiver(0:10:10*Nv, zeros(1,Nv+1),[0, squeeze(Fdyn(1,2:Nv+1,7)-Fdyn(1,1:Nv,7))],...
    [0, squeeze(Fdyn(2,2:Nv+1,7)-Fdyn(2,1:Nv,7))]);
box off;
xlim([0 200]);
ylim([-10 30]);
subplot(3,2,4);
quiver(0:10:10*Nv, zeros(1,Nv+1),[0, squeeze(Fdyn_(1,2:Nv+1,7)-Fdyn_(1,1:Nv,7))],...
    [0, squeeze(Fdyn_(2,2:Nv+1,7)-Fdyn_(2,1:Nv,7))]);
box off;
xlim([0 200]);
ylim([-10 30]);
subplot(3,2,5);
quiver(0:10:10*Nv, zeros(1,Nv+1),squeeze(mean(PV(:,1,1:Nv+1,7)))',...
    squeeze(mean(PV(:,2,1:Nv+1,7)))');
box off;
xlim([0 200]);
ylim([-10 30]);
xlabel('time [ms]');
subplot(3,2,6);
quiver(0:10:10*Nv, zeros(1,Nv+1),squeeze(mean(PV_(:,1,1:Nv+1,7)))',...
    squeeze(mean(PV_(:,2,1:Nv+1,7)))');
box off;
xlim([0 200]);
ylim([-10 30]);
xlabel('time [ms]');
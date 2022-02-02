%% Homework 2
clc, close all, clear all

%% 1.a

load('HW2.1.mat');
fs = 1/(timestep(2)-timestep(1));
T = 1/fs;

%% 1.a.3

angles = zeros(size(G));
for n = 2:length(angles)
    angles(:, n) = angles(:, n-1) + G(:, n)/fs;
end

for n = 0:0.001:3.14
    aa = [cos(n) 0 sin(n); 0 1 0 ; -sin(n) 0 cos(n)]*angles;
    avgang(round(1+1000*n)) = mean(abs(aa(1, :)));
end
[~, rotang1] = min(avgang);
rotang1 = rotang1/1000;
rotmat1 = [cos(rotang1) 0 sin(rotang1); 0 1 0 ; -sin(rotang1) 0 cos(rotang1)];
angles = rotmat1*angles;
angles = mod(angles+pi, 2*pi) - pi;
G_rot = rotmat1*G;

figure(3); hold on;
plot(timestep, angles')
xlabel('Time [s]'); ylabel('Angle [rad]')
legend('x-axis orientation', 'y-axis orientation', 'z-axis orientation')

%% 1.a.1



Acc = rotmat1*A;

A_abs = sum(Acc.^2).^0.5;
A_abshorz = sum(Acc(1:2, :).^2).^0.5; % absolute value of the acceleration on the horizontal plane
A_locmean = movmean(A_abshorz, [15 15]);
A_locvariance = movmean((A_abshorz - A_locmean).^2, [15 15]);

h = [1/2 1/2];
binomialCoeff = conv(h,h);
for n = 1:150
    binomialCoeff = conv(binomialCoeff,h);
end

[peak, peakind] = findpeaks(conv(A_locvariance, binomialCoeff));
[vall, vallind] = findpeaks(-conv(A_locvariance, binomialCoeff));
vall = -vall;

N_steps = length(vall);     % this is the estimated number of steps


figure(1);
hold on;
plot(peakind, peak, 'o')
plot(vallind, vall, 'o')
plot([zeros(1, (length(binomialCoeff)-1)/2) A_locvariance])
plot(conv(A_locvariance, binomialCoeff))
xlabel('Time instants'); ylabel('Local acceleration variance [m^2 s^-4]');
legend('Local acc. var. peaks', 'Local acc. var. valleys', 'Local acc. variance', 'Smoothed loc. acc. var.')

%% 1.a.2


Atilde = lowpass(A_abshorz, 3, fs);

aux1 = movmax(Atilde, 15);
aux2 = movmin(Atilde, 15);
for n = 1:length(vall)    
    SL(n) = (aux1(vallind(n)-(length(binomialCoeff)-1)/2) - aux2(vallind(n)-(length(binomialCoeff)-1)/2))^0.25;
end
SL = SL/mean(SL)*0.8;

figure(2)
hold on;
plot(mean(SL)*ones(size(SL)), ':')
plot(SL, 'o')
plot(SL, '-.')
ylim([0, max(SL)+0.2])
legend('Avg. step length')
xlabel('Step number'); ylabel('Step length')

%% 1.a.4

positions = zeros(2, N_steps+1);
for i = 1:N_steps
    positions(1, 1+i) = positions(1, i) + SL(i)*cos(angles(3, vallind(i)-(length(binomialCoeff)-1)/2));
    positions(2, 1+i) = positions(2, i) + SL(i)*sin(angles(3, vallind(i)-(length(binomialCoeff)-1)/2));
end
radius = sum(mean(positions').^2)^0.5;  % this is the radius of the positions 
% computing the LS circle parameters
ac=[positions(1, :)' positions(2, :)' ones(size(positions(1, :)'))]\[-(positions(1, :)'.^2+positions(2, :)'.^2)];
xc = -.5*ac(1);
yc = -.5*ac(2);
R  =  sqrt((ac(1)^2+ac(2)^2)/4-ac(3));
ac=[positions(1, :)' positions(2, :)' ones(size(positions(1, :)'))]\[-(positions(1, :)'.^2+positions(2, :)'.^2)];
xc = -.5*ac(1);
yc = -.5*ac(2);
Rc  =  sqrt((ac(1)^2+ac(2)^2)/4-ac(3));
circle = [[xc+Rc.*cos(0:0.01:3.14) xc+Rc.*cos(3.14:-0.01:0)]', [yc+Rc.*sin(0:0.01:3.14) yc-Rc.*sin(3.14:-0.01:0)]']';


% plotting the positions
figure(4); hold on;
plot(positions(1, :), positions(2, :), '.-');
plot([xc+Rc.*cos(0:0.01:3.14) xc+Rc.*cos(3.14:-0.01:0)] , [yc+Rc.*sin(0:0.01:3.14) yc-Rc.*sin(3.14:-0.01:0)], ':black')
xlabel('x-axis [m]'); ylabel('y-axis [m]')
xlim([min(positions(1, :))-1, max(positions(1, :))+1])
ylim([min(positions(2, :))-1, max(positions(2, :))+1])
legend('Step by step est. trajectory')

%% 1.a.5

T = mean(timestep(vallind-(length(binomialCoeff)-1)/2)-[0; timestep(vallind(1:end-1)-(length(binomialCoeff)-1)/2)]);
P = [1 0 T 0; ...
     0 1 0 T; ...
     0 0 1 0; ...
     0 0 0 1];
H = [1 0 0 0; 0 1 0 0];
R = 1*diag([5 5]);
Ct = 10*diag([0.5 0.5 1 1]);
Cp = 10*Ct;
theta = zeros(4, 1);

for i = 2:N_steps
    % Prediction
    theta_k = P*theta;
    Cp_k = P*Cp*P' + Ct;    % variance of the predictions
    % Measurements
    Zk = H*theta_k;
    Sk = H*Cp_k*H' + R;
    % Update
    Gk = (Cp_k*H')*pinv(Sk); % KF gain
    E_k = positions(:, i) - Zk;
    theta_k1 = theta_k + Gk*E_k;
    Cp_k1 = (eye(4) - Gk*H)*Cp_k;
    
    Pest(i, :) = theta_k1;
    
    theta = theta_k1;
    Cp = Cp_k1;
    
end

RMSEtraj = 0;
for i = 1:length(Pest)
    RMSEtraj = RMSEtraj + min((sum((Pest(i, 1:2) - circle')'.^2)).^0.5)/length(Pest);
end

figure(5); hold on;
plot(positions(1, :), positions(2, :), '.-');
plot(Pest(:, 1), Pest(:, 2))
legend('Pedestrian traj.', 'KF Tracked traj.');
xlim([min(positions(1, :))-1, max(positions(1, :))+1])
ylim([min(positions(2, :))-1, max(positions(2, :))+1])

%% 2.a

ts = 1/fs;   % sampling period

sigma_ah = mean((A_abshorz - mean(A_abshorz)).^2);
sigma_g = mean((G_rot(3, :) - mean(G_rot(3, :))).^2);

theta = [0 0 1.5 0 mean(G_rot(3,:)) mean(A_abshorz)]';   % xposition -yposition - abs.speed - 
Var = diag([sigma_ah sigma_g]);
H   = [ 0 0 0 0 0 1; ...
        0 0 0 0 1 0];
F_k = zeros(length(theta));
F_k(1,1) = 1;
F_k(2,2) = 1;
F_k(1,3) = ts*cos(theta(4));
F_k(1,4) = ts*sin(theta(4));
F_k(2,3) = -(ts*theta(3) + ts^2/2*theta(6))*sin(theta(4));
F_k(2,4) =  (ts*theta(3) + ts^2/2*theta(6))*cos(theta(4));
F_k(3,6) = ts;
F_k(4,5) = ts;

Ct = zeros(length(theta));
Ct(1,1) = ts^5/20*(Var(2,2)*(theta(3)*sin(theta(4)))^2 + Var(1,1)*(theta(3)*cos(theta(4)))^2);
Ct(1,2) = ts^5/20*sin(theta(4))*cos(theta(4))*(Var(1,1) - Var(2,2)*theta(3)^2);
Ct(1,3) = ts^4/8*Var(1,1)*cos(theta(4));
Ct(1,4) = - ts^4/8*Var(2,2)*theta(3)*sin(theta(4));
Ct(1,5) = - ts^3/6*Var(2,2)*theta(3)*sin(theta(4));
Ct(1,6) = ts^3/6*Var(1,1)*cos(theta(4));
Ct(2,1) = Ct(1,2);
Ct(2,2) = ts^5/20*(Var(2,2)*(theta(3)*cos(theta(4)))^2 + Var(1,1)*(theta(3)*sin(theta(4)))^2);
Ct(2,3) = ts^4/8*Var(1,1)*sin(theta(4));
Ct(2,4) = ts^4/8*Var(2,2)*theta(3)*cos(theta(4));
Ct(2,5) = ts^3/6*Var(2,2)*theta(3)*cos(theta(4));
Ct(2,6) = ts^3/6*Var(1,1)*sin(theta(4));
Ct(3,1) = Var(1,1)*ts^4/8*cos(theta(4));
Ct(3,2) = Var(1,1)*ts^4/8*sin(theta(4));
Ct(3,3) = Var(1,1)*ts^3/3;
Ct(3,6) = Var(1,1)*ts^2/2;
Ct(4,1) = -Var(2,2)*ts^4/8*theta(3)*sin(theta(4));
Ct(4,2) = Var(2,2)*ts^4/8*theta(3)*cos(theta(4));
Ct(4,4) = Var(2,2)*ts^3/3;
Ct(4,5) = Var(2,2)*ts^2/2;
Ct(5,1) = -Var(2,2)*ts^3/6*theta(3)*sin(theta(4));
Ct(5,2) = Var(2,2)*ts^3/6*theta(3)*cos(theta(4));
Ct(5,4) = Var(2,2)*ts^2/2;
Ct(5,5) = Var(2,2)*ts;
Ct(6,1) = Var(1,1)*ts^3/6*cos(theta(4));
Ct(6,2) = Var(1,1)*ts^3/6*sin(theta(4));
Ct(6,3) = Var(1,1)*ts^2/2;
Ct(6,6) = Var(1,1)*ts;


R = diag([sigma_ah sigma_g]);

Pest = zeros(length(timestep), 2);  %estimated positions
Pest(1,:) = theta(1:2);
Vest = zeros(1, length(timestep));  % estimated velocity
Vest(1) = theta(3);
Cp = 10*Ct;
theta_k1 = theta;

for t = 2:length(timestep)
    %Prediction
    theta_k(1) = theta(1) + (theta(3)*ts + theta(6)*ts^2/2)*cos(theta(4)); %+ v(1);
    theta_k(2) = theta(2) + (theta(3)*ts + theta(6)*ts^2/2)*sin(theta(4)); %+ v(2);
    theta_k(3) = 0.8/T + theta(6)*ts;
    theta_k(4) = wrapToPi(theta(4) + theta(5)*ts); %+ v(4);
    theta_k(5) = wrapToPi(theta(5)); %+ v(5);
    theta_k(6) = theta(6); %+ v(6);
    if (size(theta_k) == [1 6])
        theta_k = theta_k';
    end
    Cp_k = F_k*Cp*F_k' + Ct;
    %Measurement prediction
    Zk = H*theta_k;
    %Predicted noise covariance
    Sk = H*Cp_k*H' + R;
    %Update
    %Kalman Gain
    G_k = (Cp_k*H')*pinv(Sk);
    E_k = [A_abshorz(t); G_rot(3,t-1)] - Zk;
    theta_k1 = theta_k + G_k*E_k;
    Cp_k1 = (eye(6) - G_k*H) * Cp_k;
    %Actual position estimation
    Pest(t,:) = theta_k1(1:2)';
    Sest(t) = theta_k1(3);
    %Update of the new prior variables and matrices
    theta = theta_k1;
    Cp = Cp_k1;
    
    F_k = zeros(length(theta));
    F_k(1,1) = 1;
    F_k(2,2) = 1;
    F_k(1,3) = ts*cos(theta(4));
    F_k(1,4) = ts*sin(theta(4));
    F_k(2,3) = -(ts*theta(3) + ts^2/2*theta(6))*sin(theta(4));
    F_k(2,4) =  (ts*theta(3) + ts^2/2*theta(6))*cos(theta(4));
    F_k(3,6) = ts;
    F_k(4,5) = ts;
    
    Ct = zeros(length(theta));
    Ct(1,1) = ts^5/20*(Var(2,2)*(theta(3)*sin(theta(4)))^2 + Var(1,1)*(theta(3)*cos(theta(4)))^2);
    Ct(1,2) = ts^5/20*sin(theta(4))*cos(theta(4))*(Var(1,1) - Var(2,2)*theta(3)^2);
    Ct(1,3) = ts^4/8*Var(1,1)*cos(theta(4));
    Ct(1,4) = - ts^4/8*Var(2,2)*theta(3)*sin(theta(4));
    Ct(1,5) = - ts^3/6*Var(2,2)*theta(3)*sin(theta(4));
    Ct(1,6) = ts^3/6*Var(1,1)*cos(theta(4));
    Ct(2,1) = Ct(1,2);
    Ct(2,2) = ts^5/20*(Var(2,2)*(theta(3)*cos(theta(4)))^2 + Var(1,1)*(theta(3)*sin(theta(4)))^2);
    Ct(2,3) = ts^4/8*Var(1,1)*sin(theta(4));
    Ct(2,4) = ts^4/8*Var(2,2)*theta(3)*cos(theta(4));
    Ct(2,5) = ts^3/6*Var(2,2)*theta(3)*cos(theta(4));
    Ct(2,6) = ts^3/6*Var(1,1)*sin(theta(4));
    Ct(3,1) = Var(1,1)*ts^4/8*cos(theta(4));
    Ct(3,2) = Var(1,1)*ts^4/8*sin(theta(4));
    Ct(3,3) = Var(1,1)*ts^3/3;
    Ct(3,6) = Var(1,1)*ts^2/2;
    Ct(4,1) = -Var(2,2)*ts^4/8*theta(3)*sin(theta(4));
    Ct(4,2) = Var(2,2)*ts^4/8*theta(3)*cos(theta(4));
    Ct(4,4) = Var(2,2)*ts^3/3;
    Ct(4,5) = Var(2,2)*ts^2/2;
    Ct(5,1) = -Var(2,2)*ts^3/6*theta(3)*sin(theta(4));
    Ct(5,2) = Var(2,2)*ts^3/6*theta(3)*cos(theta(4));
    Ct(5,4) = Var(2,2)*ts^2/2;
    Ct(5,5) = Var(2,2)*ts;
    Ct(6,1) = Var(1,1)*ts^3/6*cos(theta(4));
    Ct(6,2) = Var(1,1)*ts^3/6*sin(theta(4));
    Ct(6,3) = Var(1,1)*ts^2/2;
    Ct(6,6) = Var(1,1)*ts;

end

avgspeedEKF = mean(((Pest(2:end, 1) - Pest(1:end-1, 1)).^2 + (Pest(2:end, 2) - Pest(1:end-1, 2)).^2).^0.5/ts);

% computing the LS circle parameters
ac=[Pest(:, 1) Pest(:, 2) ones(size(Pest(:, 1)))]\[-(Pest(:, 1).^2+Pest(:, 2).^2)];
xc = -.5*ac(1);
yc = -.5*ac(2);
Rc  =  sqrt((ac(1)^2+ac(2)^2)/4-ac(3));
circle = [[xc+Rc.*cos(0:0.01:3.14) xc+Rc.*cos(3.14:-0.01:0)]', [yc+Rc.*sin(0:0.01:3.14) yc-Rc.*sin(3.14:-0.01:0)]']';

RMSEekf = 0;
for i = 1:length(Pest)
    RMSEekf = RMSEekf + min((sum((Pest(i, :) - circle')'.^2)).^0.5)/length(Pest);
end

figure(6); hold on;
plot(positions(1, :), positions(2, :), '.-');
plot(Pest(:, 1), Pest(:, 2))
plot([xc+Rc.*cos(0:0.01:3.14) xc+Rc.*cos(3.14:-0.01:0)] , [yc+Rc.*sin(0:0.01:3.14) yc-Rc.*sin(3.14:-0.01:0)], ':black')
legend('Pedestrian traj.', 'EKF Tracked traj.');
xlim([min(positions(1, :))-2, max(positions(1, :))+1])
ylim([min(positions(2, :))-1, max(positions(2, :))+2])

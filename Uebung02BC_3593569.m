clear all
close all
clc
% Uebung 2  Aufgabe a  BC
% Erdrotation
w_E = [0;0;7.292115e-5];
Omega_iee = [0 -w_E(3) 0
             w_E(3) 0 0
             0 0 0];

% 读取IMU数据 ax ay az(m/s^2) gx gy gz(rad/s)
IMU.raw = csvread('uebung01-IMU.csv', 1, 0);
t = IMU.raw(:,1);
IMUdata(:,1:3) = IMU.raw(:,2:4);
IMUdata(:,4:6) = IMU.raw(:,5:7);
REF.raw = csvread('uebung01-REF.csv', 1, 0);
refp = REF.raw(:,2:4);
refv = REF.raw(:,5:7)';

%角度转弧度 deg2rad
refRPY = deg2rad(REF.raw(:,8:10));
refp(:,1:2) = deg2rad(refp(:,1:2)); 

%Ref Position in e-system
wgs84 = wgs84Ellipsoid('meter');
[xRef(:),yRef(:),zRef(:)] = geodetic2ecef(refp(:,1),refp(:,2),refp(:,3),wgs84);

%DCM t0 berechnen
Cne_0 = C(3,-refp(1,2))*C(2,refp(1,1)+pi/2);
Cbn_0 = C(3,-refRPY(1,3))*C(2,-refRPY(1,2))*C(1,-refRPY(1,1));
Cpe_0 = Cne_0*Cbn_0;
Omega_iep = inv(Cpe_0)*Omega_iee*Cpe_0;
w_ieP(1,1) = Omega_iep(3,2);
w_ieP(2,1) = Omega_iep(1,3);
w_ieP(3,1) = Omega_iep(2,1);
refv(:,1) = Cne_0*refv(:,1);           % v in t0                                      

%DCM t1 berechnen
Cne_1 = C(3,-refp(2,2))*C(2,refp(2,1)+pi/2);
Cbn_1 = C(3,-refRPY(2,3))*C(2,-refRPY(2,2))*C(1,-refRPY(2,1));
Cpe_1 = Cne_1*Cbn_1;
Omega_iep = inv(Cpe_1)*Omega_iee*Cpe_1;
w_ieP(1,1) = Omega_iep(3,2);
w_ieP(2,1) = Omega_iep(1,3);
w_ieP(3,1) = Omega_iep(2,1);
refv(:,2) = Cne_1*refv(:,2);           % v in t1                                       

% ref pos 
x(:,1) = xRef(:);
x(:,2) = yRef(:);
x(:,3) = zRef(:);
                                                

%Quaternion t0 im e-system berechnen 
qt0 = compact(quaternion((Cpe_0),'rotmat','frame'));

%Quaternion t1 im e-system berechnen
qt1 = compact(quaternion((Cpe_1),'rotmat','frame'));

%DCM Brechnen

for i = 3:length(t)
    Cpe{1} = Cpe_0;
    Cpe{2} = Cpe_1;
    q(1,:) = qt0;
    q(2,:) = qt1;
    delta_t = abs(t(i-1) - t(i-2));

    %3.4
    delta_alpha(i-1,1:3) = (IMUdata(i-1,4:6)+IMUdata(i-2,4:6))*delta_t/2;
    delta_alpha(i,1:3) = (IMUdata(i,4:6)+IMUdata(i-1,4:6))*delta_t/2;

    %3.9
    delta_beta(i-1,:) = delta_alpha(i-1,1:3)'-Cpe{i-2}'*w_E*delta_t;
    delta_beta(i,:) = delta_alpha(i,1:3)'-Cpe{i-1}'*w_E*delta_t;

    %3.15
    omegaepp(i-2,:) = (3*delta_beta(i-1,:)-delta_beta(i,:))/(2*delta_t);
    omegaepp(i-1,:) = (delta_beta(i-1,:)+delta_beta(i,:))/(2*delta_t);
    omegaepp(i,:) = (3*delta_beta(i,:)+delta_beta(i-1,:))/(2*delta_t);


    %3.2 RK3
    [Cpe{i},q(i,:)] = RK3Quaterionen(@TimeDerivativeQuaterionen,omegaepp(i-2:i,:),q(i-2,:),2*delta_t);

    %4.3 ΔV sowie 3.4
    delta_v(i-1,:) =  (IMUdata(i-1,1:3)+IMUdata(i-2,1:3))*delta_t/2;
    delta_v(i,:) = (IMUdata(i,1:3)+IMUdata(i-1,1:3))*delta_t/2;

    %4.12 V
    v_e(1,:) = refv(:,1);
    v_e(2,:) = refv(:,2);
     if i >= 3 
             v_e(i,:) = v_e(i-2,:)'+(Cpe{i-2}*(3*delta_v(i-1,:)-delta_v(i,:))')+ ... 
                          4*Cpe{i-1}*(delta_v(i-1,:)' + delta_v(i,:)') + ...
                          Cpe{i}*(3*delta_v(i,:)' - delta_v(i-1,:)')/6 - ...
                          (2*Omega_iee*refv(:,i-2)+ ...
                          Omega_iee*Omega_iee*x(i-2,:)')*2*delta_t;
             Aufgabe = 1;%控制最后图片的名字 Control pic name
             %Aufgabe 2 
                  if i == 8
                       V_fehler = [0; 2; 0];
                       v_e(i,:) = v_e(i,:) + V_fehler'; 

                  end
             Aufgabe = 2;
     end

    %4.15 x
    x_v0304(1,:) = x(1,:);
    x_v0304(2,:) = x(2,:);
    x_v0304(i,:) = x(i-1,:) + v_e(i-1,:)*delta_t;

end




Rk3_temp = [x(1,:),refv(:,1)',q(1,:)];
for i = 2:length(t)
    x0 = [IMUdata(i,1:3) IMUdata(i,4:6)];
    Rk3_temp(i,:) = RungeKutta3(@TimeDerivativePosVelAtt_e,Rk3_temp(i-1,:)',x0',delta_t);

    %Aufgabe b 
        if i == 5
                    error_vector = [0; 2; 0; 0; 0; 0; 0; 0; 0; 0]; 
                    Rk3_temp(i,1:10) = Rk3_temp(i,1:10) + error_vector';

        end
         
end
RK3_Pos = Rk3_temp(:,1:3);
RK3_Vel = Rk3_temp(:,4:6);
RK3_q = Rk3_temp(:,7:10);

for i = 1:length(RK3_q)
    RK3_rotm{i,1} = rotmat(quaternion(RK3_q(i,:)),'frame');
end







%% plot
figure
plot3(x(:,1),x(:,2),x(:,3))
hold on
plot3(x_v0304(:,1),x_v0304(:,2),x_v0304(:,3))
plot3(RK3_Pos(:,1),RK3_Pos(:,2),RK3_Pos(:,3))
xlabel("[m]");ylabel("[m]");zlabel("[m]");
title(['Aufgabe ', num2str(Aufgabe), ' 3D Plot for all']);
legend('REF', 'V0304', 'RK3');
saveas(gcf, ['Aufgabe ', num2str(Aufgabe), ' 3D Plot for all.png']);

figure
plot(x(:,2),x(:,3));
hold on 
plot(x_v0304(:,2),x_v0304(:,3));
plot(RK3_Pos(:,2),RK3_Pos(:,3));
title(['Aufgabe ', num2str(Aufgabe), '.1 YZ-ECEF Plot for all']);
xlabel("[m]");ylabel("[m]");
legend('REF', 'V0304', 'RK3');
saveas(gcf, ['Aufgabe ', num2str(Aufgabe), '.1 YZ-ECEF Plot for all.png']);


%%
[llax] = ecef2lla(x);
[llax0304] = ecef2lla(x_v0304);
[llaxrk3] = ecef2lla(RK3_Pos);
figure
plot3(llax(:,1),llax(:,2),llax(:,3))
hold on
plot3(llax0304(:,1),llax0304(:,2),llax0304(:,3))
plot3(llaxrk3(:,1),llaxrk3(:,2),llaxrk3(:,3))
xlabel("[m]");ylabel("[m]");zlabel("[m]");
title(['Aufgabe ', num2str(Aufgabe), '.2 3D Plot for all']);
legend('REF', 'V0304', 'RK3');
saveas(gcf, ['Aufgabe ', num2str(Aufgabe), '.2 3D Plot for all.png']);

figure
plot(llax(:,3));
hold on
plot(llax0304(:,3));
plot(llaxrk3(:,3));
ylabel("[m]")
title(['Aufgabe ', num2str(Aufgabe), '.2 Ellipsoidheight Plot for all']);
legend('REF', 'V0304', 'RK3');
saveas(gcf, ['Aufgabe ', num2str(Aufgabe), '.2 Ellipsoidheight Plot for all.png']);
%%

for i = 1:length(q)
eulerWinkel = euler(quaternion(q(i,:)),"YZX","frame");
yaw0304(i) = rad2deg(eulerWinkel(1));
pitch0304(i) = rad2deg(eulerWinkel(2));
roll0304(i) = rad2deg(eulerWinkel(3));
end 
for i = 1:length(Rk3_temp(:,7:10))
eulerWinkel = euler(quaternion(Rk3_temp(i,7:10)),"YZX","frame");
yawrk3(i) = rad2deg(eulerWinkel(1));
pitchrk3(i) = rad2deg(eulerWinkel(2));
rollrk3(i) = rad2deg(eulerWinkel(3));
end 
refRPY = REF.raw(:,8:10);

figure
subplot(3,1,1); 
plot(refRPY(:,1));
hold on
plot(roll0304);
plot(rollrk3);
ylabel("[deg]")
legend('REF', 'V0304', 'RK3');
title('ROLL');

subplot(3,1,2); 
plot(refRPY(:,2));
hold on
plot(pitch0304);
plot(pitchrk3);
ylabel("[deg]")
legend('REF', 'V0304', 'RK3');
title('PITCH');

subplot(3,1,3);
plot(refRPY(:,3));
hold on
plot(yaw0304);
plot(yawrk3);
ylabel("[deg]")
legend('REF', 'V0304', 'RK3');
title('YAW');
sgtitle(['Aufgabe ', num2str(Aufgabe), '.3 RPY Plot for all']);
saveas(gcf, ['Aufgabe ', num2str(Aufgabe), '.3 RPY Plot for all.png']);
%%
refv = REF.raw(:,5:7);
for i = 1:length(refp)
Cne = C(3,-refp(i,2))*C(2,refp(i,1)+pi/2);
Cen = Cne';
v_0304_ned(i,:) = (Cen*v_e(i,:)')';
v_rk3_ned(i,:) = (Cen*RK3_Vel(i,:)')';
end
figure
subplot(3,1,1); 
plot(refv(:,1),'o');
hold on
plot(v_0304_ned(:,1));
plot(v_rk3_ned(:,1));
ylabel("[m/s]")
legend('REF', 'V0304', 'RK3');
title('Nord');

subplot(3,1,2); 
plot(refv(:,2),'o');
hold on
plot(v_0304_ned(:,2));
plot(v_rk3_ned(:,2));
ylabel("[m/s]")
legend('REF', 'V0304', 'RK3');
title('Ost');

subplot(3,1,3);
plot(refv(:,3),'o');
hold on
plot(v_0304_ned(:,3));
plot(v_rk3_ned(:,3));
ylabel("[m/s]")
legend('REF', 'V0304', 'RK3');
title('Down');
sgtitle(['Aufgabe ', num2str(Aufgabe), '.4 NED Geschwindigkeit Plot for all']);
saveas(gcf, ['Aufgabe ', num2str(Aufgabe), '.4 NED Geschwindigkeit Plot for all.png']);
%% Funktion
function y = RungeKutta3(func,y0,x0,t)
k1 = func(y0,x0);
k2 = func(y0+t/2*k1,x0);
k3 = func(y0-t*k1+2*t*k2,x0);
y = y0 + t/6*(k1+4*k2+k3);
y(7:10) = compact(normalize(quaternion(y(7:10)')));
end

function [rotm,quat] = RK3Quaterionen(func,omega,q0,h)

k1 = func(omega(1,:),q0);
k2 = func(omega(2,:),q0+k1*h/2);
k3 = func(omega(3,:),q0-k1*h+k2*2*h);
q = normalize(quaternion(q0 + h/6*(k1+4*k2+k3)));
quat = compact(q);
rotm = rotmat(q,'frame');
end

function qdot = TimeDerivativeQuaterionen(omega_epp,q_pe)
w1 = omega_epp(1);
w2 = omega_epp(2);
w3 = omega_epp(3);

A = [0 w1 w2 w3
    -w1 0 w3 -w2
    -w2 -w3 0 w1
    -w3 w2 -w1 0];
qdot = (1/2*A*q_pe')';
end

function [y_dot] = TimeDerivativePosVelAtt_e(y,x)
   w_E = [0;0;7.292115e-5];
   Omega_iee = [0 -w_E(3) 0
             w_E(3) 0 0
             0 0 0];
 

  % Data extraction
  pos_e = y(1:3);
  v_e = y(4:6);
  
  C_p2e_quat = y(7:10);
  C_p2e_rotm = rotmat(quaternion(C_p2e_quat'),'frame');

  a_ip_p = x(1:3);
  omega_ip_p = x(4:6);

  % Calculations
  omega_ep_p = omega_ip_p - C_p2e_rotm' * w_E;
  A = [      0.0     , omega_ep_p(1), omega_ep_p(2), omega_ep_p(3);...
        -omega_ep_p(1),       0.0     ,  omega_ep_p(3), -omega_ep_p(2);...
        -omega_ep_p(2), -omega_ep_p(3),       0.0     ,  omega_ep_p(1);...
        -omega_ep_p(3),  omega_ep_p(2), -omega_ep_p(1),       0.0     ];

  y_dot = [v_e
           C_p2e_rotm * a_ip_p - 2 * Omega_iee * v_e -...
           Omega_iee * Omega_iee * pos_e
           0.5 * A * C_p2e_quat];
end




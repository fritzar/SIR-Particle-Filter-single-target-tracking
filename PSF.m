clear
clc
close all
format long;
L = 0; %？
Monte = 50;
axisX = 30; %50%x轴最大值
axisY = 30;%50
Total_time = 20;%40
Re_x = 1; %Resolution of x-axis
Re_y = 1;
ypart = 0:Re_y:axisY; %步长为分辨率
xpart = 0:Re_x:axisX;
numY = length(ypart); %取长度
numX = length(xpart);
Sigma_noise = 1; %?
SNR = [6,9]; %
NpN =100;%2.^[1:10]--仿真粒子数变化
%Np = 1000;
T_step = 1; % The size of the time cell:Time_step
q1 = 0.0015; %"空间过程噪声加速度的功率谱密度 ?

F = [1 T_step 0 0 %状态转移矩阵
    0   1    0 0
    0   0    1 T_step
    0   0    0 1];
Q=[T_step^3*q1/3  T_step^2*q1/2  0           0 ;
    T_step^2*q1/2 T_step*q1      0           0 ;
    0             0           T_step^3*q1/3  T_step^2*q1/2;
    0             0           T_step^2*q1/2  q1*T_step];        % ProcessNoise covariance matrix

Target_number = 1; %单目标
velocity_init = 1; %
[initx,x,x_c] = GenerateTarget(Target_number,velocity_init,axisX,axisY,Total_time,F,Q);
% x = [linspace(10,25,Total_time);linspace(9,21 ,Total_time)];
x_dis = ceil(x(1,:)/Re_x)*Re_x; %能分辨的目标位置， ceil朝正无穷方向取整
y_dis = ceil(x(3,:)/Re_y)*Re_y;

figure(1)
plot(x(1,:),x(3,:),'bp-');
axis equal
axis([0 30 0 30])
hold on; grid on;

E_target = zeros(7,Total_time,Monte,length(SNR),length(NpN));
Target_p_error = zeros(Total_time,Monte,length(SNR),length(NpN));
error_P = zeros(Total_time,length(SNR),length(NpN));
for Np_i = 1:length(NpN)
    Np = NpN(Np_i);
    for snr_i = 1:length(SNR)
        snr = SNR(snr_i);
        % Frame_data = zeros(numY,numX,Total_time);
          %% 高斯噪声 %%%
            I_noise=Sigma_noise*randn(numY,numX,Total_time);
            Q_noise=Sigma_noise*randn(numY,numX,Total_time);
            %
            Frame_data=(1/2^0.5).*(I_noise+Q_noise.*(-1)^0.5);
         %% 瑞利噪声
            % Noise_data = raylrnd(Sigma_noise,numY,numX*Total_time);
        for monte_i = 1:Monte
            display(['Np=',num2str(Np),'; Monte=',num2str(monte_i)])
            %%
            xy_data = zeros(size(Frame_data));
            A=sqrt(10.^(snr/10)*Sigma_noise);
            for t = 1:Total_time
                B=A*exp(sqrt(-1)*2*pi*rand);
                %     Frame_data(:,:,t) = Noise_data(:,numX*(t-1)+1:numX*t);
                %     for m = 1:numY
                %         for n=1:numX
                %             Frame_data(m,n,t) = Frame_data(m,n,t) + B*exp(-(((n*Re_x-x_dis(t))^2/(2*Re_x)+(m*Re_y-y_dis(t))^2/(2*Re_y))*L));
                %         end
                %     end
                %     Frame_data(y_dis(t),x_dis(t),t) = raylrnd(sqrt(Sigma_noise+A^2),1);
                Frame_data(y_dis(t),x_dis(t),t) = Frame_data(y_dis(t),x_dis(t),t)+ B;
                xy_data(:,:,t)=abs(Frame_data(:,:,t));
                
                %     figure(1)
                %     imagesc(xpart,ypart,xy_data(:,:,t));
                %     colorbar
                %     set(gca,'yDir','normal');
                %     axis equal
                %     axis tight
                % %     L = L/2;
                % %     surf(xy_data(:,:,t));
                %     pause(0.5)
                %     gif_f=getframe(gcf);
                % imind=frame2im(gif_f);
                % [imind,cm] = rgb2ind(imind,256);
                % if t==1
                %     imwrite(imind,cm,'non-psf plane','gif', 'Loopcount',inf,'DelayTime',0.5);%第一次必须创建！
                % else
                %     imwrite(imind,cm,'non-psf plane','gif','WriteMode','append','DelayTime',0.5);
                % end
                %     pause
            end
            
            [E_target_state]=standard_SIR(Np,initx,Re_x,Re_y,numX,numY,Total_time,xy_data,Sigma_noise,A);
            E_target(:,:,monte_i,snr_i,Np_i) = E_target_state;
        end
Target_p_error(:,:,snr_i,Np_i) = (squeeze(E_target(1,:,:,snr_i,Np_i))-repmat(x(1,:),Monte,1)').^2 + (squeeze(E_target(4,:,:,snr_i,Np_i))-repmat(x(3,:),Monte,1)').^2; %T*Np*length(SNR)
error_P(:,snr_i,Np_i) = squeeze(sqrt(mean(Target_p_error(:,:,snr_i,Np_i),2))); %%T*length(SNR)
  end

xy_P = E_target(:,:,ceil(Monte*rand),1,Np_i);
end;
%xy_P3 = E_target(:,:,ceil(Monte*rand),snr_i,Np_i); %最后一次（snr_i和Np_i停在最后一个值）ceil(Monte*rand)

figure(2)
plot(xy_P(1,:),xy_P(4,:),'bp-')
axis([0,axisX,0,axisY])
hold on;grid on;
plot(x(1,:),x(3,:),'ko-')
title('跟踪结果')
xlabel('x方向距离')
ylabel('y方向距离')
legend('估计轨迹','真实轨迹')

% figure(3)
% plot(error_P(:,1,Np_i),'^-');
% title('各帧均方误差')
% % axis([0,Total_time,0,1])
% xlabel('时间/帧')
% ylabel('均方误差')
% grid on
% 
% figure(4)
% plot(xy_P3(1,:),xy_P3(4,:),'bp-')
% axis([0,axisX,0,axisY])
% hold on;grid on;
% plot(x(1,:),x(3,:),'ko-')
% title('跟踪结果')
% xlabel('x方向距离')
% ylabel('y方向距离')
% legend('估计轨迹','真实轨迹')
% 
% figure(5)
% plot(error_P(:,3,Np_i),'^-');
% title('各帧均方误差')
% % axis([0,Total_time,0,1])
% xlabel('时间/帧')
% ylabel('均方误差')
% grid on
% 
% figure(10)
% plot(SNR,mean(error_P,1),'ko-')
% grid on
% title('RMSE随信噪比曲线')
% xlabel('SNR/dB')
% ylabel('RMSE')
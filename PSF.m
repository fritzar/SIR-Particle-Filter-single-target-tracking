clear
clc
close all
format long;
L = 0; %？
Monte = 50;
axisX = 30; %x轴最大值
axisY = 30;
Total_time = 20;
Re_x = 1; %Resolution of x-axis
Re_y = 1;
ypart = 0:Re_y:axisY; %步长为分辨率
xpart = 0:Re_x:axisX;
numY = length(ypart); %取长度
numX = length(xpart);
Sigma_noise = 1; %?
SNR = [6]; %
NpN = [64,256,512];%2.^[1:10]--仿真粒子数变化
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


%% adding fake target 
% [fake_initx,fake_x,fake_x_c] = GenerateFakeTarget(Target_number,velocity_init,axisX,axisY,Total_time,F,Q);
% 
% figure(1)
% plot(x(1,:),x(3,:),'bp-','Linewidth',2);
% hold on; grid on;
% plot(fake_x(1,[5:15]),fake_x(3,5:15),'m-o','Linewidth',2);
% axis equal
% axis ([0 30 0 30])
% hold off
% 
% fake_x_dis = ceil(fake_x(1,[5:15])/Re_x)*Re_x; %能分辨的目标位置， ceil朝正无穷方向取整
% fake_y_dis = ceil(fake_x(3,[5:15])/Re_y)*Re_y;
%%


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
            
%             I_noise=Sigma_noise*randn(numX,numY,Total_time);
%             Q_noise=Sigma_noise*randn(numX,numY,Total_time);
            Frame_data=(1/2^0.5).*(I_noise+Q_noise.*(-1)^0.5);
            
             
%             figure(14)  
%             imagesc(ypart,xpart,abs(Frame_data(:,:,1)));
%             axis xy;
%             colormap(gray(64))
%             colorbar('YTickLabel',{' '});
%             caxis([0 3.5])
%             axis equal
%             axis tight
%             hold off
            
        %% 瑞利噪声
%             Noise_data = raylrnd(Sigma_noise,numY,numX*Total_time);
%             for ti=1:Total_time
%                 for xi=1:numX
%                     Frame_data (:,xi,ti) = Noise_data (:,xi*ti); 
%                 end
%             end
%             save RayNoise30_6dB Frame_data
            
%             figure(14)  
%             imagesc(ypart,xpart,abs(Frame_data(:,:,1)));
%             axis xy;
%             colormap(gray(64))
%             colorbar('YTickLabel',{' '});
%             caxis([0 3.5])
%             axis equal
%             axis tight
%             hold off
                     
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
                
                %% fake data adding to the plane
%                 if t>=5 && t <=15
%                     Frame_data(fake_y_dis(t-4),fake_x_dis(t-4),t) = Frame_data(fake_y_dis(t-4),fake_x_dis(t-4),t)+ B;
%                 end
                %%
                
                xy_data(:,:,t)=abs(Frame_data(:,:,t));
                
                
%                 figure(15)
%                 imagesc(ypart,xpart,xy_data(:,:,1));
%                 axis xy;
%                 colormap(gray(64))
%                 colorbar('YTickLabel',{' '})
%                 caxis([0 3.5])
%                 axis equal
%                 axis tight
%                 hold off
                
                
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
            
            [E_target_state]=standard_SIR(Np,initx,Re_x,Re_y,numX,numY,Total_time,xy_data,Sigma_noise,A,F,Q);
            E_target(:,:,monte_i,snr_i,Np_i) = E_target_state;
        end
Target_p_error(:,:,snr_i,Np_i) = (squeeze(E_target(1,:,:,snr_i,Np_i))-repmat(x(1,:),Monte,1)').^2 + (squeeze(E_target(4,:,:,snr_i,Np_i))-repmat(x(3,:),Monte,1)').^2; %T*Np*length(SNR)
error_P(:,snr_i,Np_i) = squeeze(sqrt(mean(Target_p_error(:,:,snr_i,Np_i),2))); %%T*length(SNR)*length（Np）
    end
xy_P = E_target(:,:,ceil(Monte*rand),1,Np_i); %随机选取最后一个Np_i的一次蒙特卡洛仿真
end

%% 单目标一次monte轨迹（小图）%cell模式 ctrl+enter执行，相当于命令窗口
figure(50)
plot(xy_P(1,:),xy_P(4,:),'bp-')
axis([0,axisX,0,axisY])
hold on;grid on;
plot(x(1,:),x(3,:),'ko-')
title('跟踪结果')
xlabel('x方向距离')
ylabel('y方向距离')
legend('估计点迹','真实点迹')
%小图参数
% set(gcf,'Position',[100 100 260 220]);
% set(gca,'Position',[.13 .17 .80 .74]);
% figure_FontSize=8;
% set(get(gca,'XLabel'),'FontSize',figure_FontSize,'Vertical','top');
% set(get(gca,'YLabel'),'FontSize',figure_FontSize,'Vertical','middle');
% set(findobj('FontSize',10),'FontSize',figure_FontSize);
% set(findobj(get(gca,'Children'),'LineWidth',0.5),'LineWidth',2);
%% 单目标一次monte轨迹（小图）%cell模式 ctrl+enter执行，相当于命令窗口
figure(51)
xy_P1 =  E_target(:,:,ceil(Monte*rand),1,1);
plot(xy_P1(1,:),xy_P1(4,:),'bp-')
axis([0,axisX,0,axisY])
hold on;grid on;
plot(x(1,:),x(3,:),'ko-')
title('跟踪结果')
xlabel('x方向距离')
ylabel('y方向距离')
legend('估计点迹','真实点迹')
%% 单目标多Np_i绘轨迹
% figure(51)
% plot(xy_P(1,:),xy_P(4,:),'gp-')
% plot(xy_P(1,:),xy_P(4,:),'kp-')
% plot(xy_P(1,:),xy_P(4,:),'bp-')
% axis([0,axisX,0,axisY])
% hold on;grid on;
% plot(x(1,:),x(3,:),'ko-')
% title('跟踪结果')
% xlabel('x方向距离')
% ylabel('y方向距离')
% legend('估计点迹','真实点迹')
%% 多Np_i/SNR_i绘RMSE
colorParticle={'b.','y.','g.','k.';'g^-','k^-','b^-','y^-';'bo','ro','mo','go'};
figure(60)
for Np_i = 1: length(NpN)
plot(error_P(:,1,Np_i),colorParticle{2,Np_i});
hold on
end
% for SNR_i = 1: length(SNR)
% plot(error_P(:,SNR_i,1),colorParticle{2,SNR_i});
% hold on
% end
hold off
title('各帧均方误差')
% axis([0,Total_time,0,1])
xlabel('时间/帧')
ylabel('均方误差')
grid on
%% 
figure(5)
plot(error_P(:,1,Np_i),'^-'); %某一个snr条件下的rmse
title('各帧均方误差')
% axis([0,Total_time,0,1])
xlabel('时间/帧')
ylabel('均方误差')
grid on
%小图参数
% set(gcf,'Position',[100 100 260 220]);
% set(gca,'Position',[.13 .17 .80 .74]);
% figure_FontSize=8;
% set(get(gca,'XLabel'),'FontSize',figure_FontSize,'Vertical','top');
% set(get(gca,'YLabel'),'FontSize',figure_FontSize,'Vertical','middle');
% set(findobj('FontSize',10),'FontSize',figure_FontSize);
% set(findobj(get(gca,'Children'),'LineWidth',0.5),'LineWidth',2);
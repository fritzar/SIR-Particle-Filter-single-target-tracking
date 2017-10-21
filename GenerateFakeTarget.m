function [initx,x,x_c] = GenerateFakeTarget(Target_number,velocity_init,Num_Cell_x,Num_Cell_y,Total_time_data,F,Q)
%生成目标
% initx=zeros(4,Target_number); %[0;0;0;0;] ?why 4?
% x = zeros(4,Target_number*Total_time_data); %每Total_time_data个列存放一个目标的各帧状态
% x_c = zeros(3*Target_number,Total_time_data); %why 3?
% F = [1 T_step 0 0
%      0   1    0 0
%      0   0    1 T_step
%      0   0    0 1];
% Q=[T_step^3*q1/3  T_step^2*q1/2  0           0 ;
%    T_step^2*q1/2 T_step*q1      0           0 ;
%    0             0           T_step^3*q1/3  T_step^2*q1/2;
%    0             0           T_step^2*q1/2  q1*T_step];        % ProcessNoise covariance matrix
Sigma_noise=1; %sigma_noise?
i_p = [14,2]';
%i_p = [10,6;6,Num_Cell_y/2;Num_Cell_x/2,Num_Cell_y-6;Num_Cell_x-6,10]'; %position i_p(4*2)?
v1=(0.5+0.5*rand(1,4))*velocity_init;  %1*4Columns
v = repmat(v1,2,1); %replicate v1 by 2*1
alpha = [pi/4-0.2,0.1,-pi/2-0.1,3*pi/4+0.3];
i_v = v.*[cos(alpha(1)),sin(alpha(1));cos(alpha(2)),sin(alpha(2));cos(alpha(3)),sin(alpha(3));cos(alpha(4)),sin(alpha(4))]'; %velocity
initx([1,3],:) = i_p(:,1:Target_number);
initx([2,4],:) = i_v(:,1:Target_number);

processNoise =Sigma_noise*sample_gaussian(zeros(length(Q),1),Q,11)';% Note ' here! state_noise_samples is: ss-Total_time_allfor t=1:Total_time_all
for t = 5:15
    for n = 1:Target_number
        if t==5
            x(:,t+(n-1)*11) = initx(:,n) + processNoise(t);
            x_c([1,2]+(n-1)*3,t) = x([1,3],t+(n-1)*11);
        else
            x(:,t+(n-1)*11) = F*x(:,t-1+(n-1)*11) + processNoise(t);
            x_c([1,2]+(n-1)*3,t) = x([1,3],t+(n-1)*11);
        end
    end
end




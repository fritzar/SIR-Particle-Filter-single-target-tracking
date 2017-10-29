function [ Particle_T ] = sample_generic( Particle_O,T_step,q1)
%SAMPLE_RB Summary of this function goes here
%   Detailed explanation goes here
        mean_temp_P=[Particle_O(1)+T_step*Particle_O(2);Particle_O(4)+T_step*Particle_O(5)];
        covar_temp_P=(T_step.^2).*[Particle_O(3),0;0, Particle_O(6)]+ [T_step^3*q1/3,0; 0,T_step^3*q1/3;];
        particle_temp=sample_gaussian(mean_temp_P,covar_temp_P,1);
        Particle_T(1)=particle_temp(1);
        Particle_T(4)=particle_temp(2);

%% -- 单独的速度更新
% mean_temp_V = ;
% covar_temp_V =
%  
%     

%% -- 如何直接带入状态方程计算（如同tutorial配套代码例程的sys）
       

end
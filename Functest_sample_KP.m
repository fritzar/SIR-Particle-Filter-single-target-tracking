%test function sample_KP
clear
clc
T_step=1;
q1 = 0.0015;
F = [1 T_step 0 0 %×´Ì¬×ªÒÆ¾ØÕó
    0   1    0 0
    0   0    1 T_step
    0   0    0 1];
Q=[T_step^3*q1/3  T_step^2*q1/2  0           0 ;
    T_step^2*q1/2 T_step*q1      0           0 ;
    0             0           T_step^3*q1/3  T_step^2*q1/2;
    0             0           T_step^2*q1/2  q1*T_step];        % ProcessNoise covariance matrix
Particle = [3,1,2,2,0.1,0.1];

Particle_next  = sample_KP( Particle,F,Q);

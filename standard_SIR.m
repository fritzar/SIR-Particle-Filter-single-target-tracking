function [E_target_state_MC]=standard_SIR(Np,initx,Re_x,Re_y,numX,numY,Total_time,xy_data,Sigma_noise,A)

format long;

%% ======The initialization of the basic parameters=====

T_step=1;          % The size of the time cell:Time_step
q1=0.002;          % q1,q2 the level of process noise in target motion

%% ---------- initial distribution of target state
delta_p=0.1;  % the stardard deviation of the inital distribution of position
new_velocity_Variance=0.01;             % standard deviation 0.1;
% q2=2;                  % q2 the measurement noise

%% ---------- Rao-balckwellision parameters
Q_l=q1*[T_step,0;0,T_step];
Q_n=q1*[T_step^3/3,0;0,T_step^3/3];
Q_ln=q1*[T_step^2/2,0;0,T_step^2/2];
A_1_t=eye(2)-Q_ln/(Q_n)*T_step;
Q_1_l=Q_l-(Q_ln.')*(inv(Q_n))*Q_ln;

% Target_number = 1;
E_target_state_MC=zeros(7,Total_time);
% single_run_time=zeros(Total_time);
% figure(12);hold on;plot(squeeze(Pre_T_particle(1,t,:,:)),squeeze(Pre_T_particle(4,t,:,:)),'c.',x_c(1:3:13,t),x_c(2:3:14,t),'rx','LineWidth',2,'MarkerSize',8);

%% ===================== Particle filtering =========================

%% ============== PF implementation ================
Pre_T_particle=zeros(7,Total_time,Np);            % Pre_track particles
% Pre_T_particle_ori=zeros(7,Total_time,Np);      % Pre_track particles
Pre_T_life_index=zeros(Np,1);                   % life length of the pre-tracks
Pre_T_life_quality=zeros(Np,1);                 % quality of life length of the pre-tracks
% Pre_weight0=zeros(Np,Total_time);               % Particles weights of Pre-track PF
% Pre_weight=zeros(Np,Total_time);                % Normalized Particles weights
% Pre_w_likehood_all=zeros(Np,Total_time);        % likehood part of the pre-PF weights
% Pre_track_bias=zeros(Np,Total_time);            % weight bias part of the pre-PF weights


for t = 1:Total_time
%     display(['Np=',num2str(Np),'; t=',num2str(t)]);
    %singerun_start=tic;
    %% --------------- detection procedure ----------------
    Detection_frame=xy_data(:,:,t);
%     clean_Frame = zeros(size(Detection_frame));
    if t==1
        %disp(['t==1']);
        % number of the birth pre-tracks
        index_x=initx(1)/Re_x;
        index_y=initx(3)/Re_y;
        index_vx=initx(2)/Re_x;
        index_vy=initx(4)/Re_y;
        % -------generate the new partitions of particles
        %--------generate position based the detection measurements
        position_x_p=repmat(index_x,Np,1)+delta_p*randn(Np,1);
        position_y_p=repmat(index_y,Np,1)+delta_p*randn(Np,1);
        %% --------generate velocity based on the detections
        velocity_x_p=repmat(index_vx,Np,1);
        velocity_y_p=repmat(index_vy,Np,1);
%        velocity_x_p=ones(Np,1);
%        velocity_y_p=ones(Np,1);
        %--------generate velocity variance
        velocity_p_kk1=new_velocity_Variance.*ones(Np,1);
        
        %--------new_pretrack=zeros(4,Num_n_target,Np);
        Pre_T_life_index=ones(1,Np);
        Pre_T_life_quality=ones(1,Np);
        for i=1:Np
            Pre_T_particle(1:6,t,i)=[position_x_p(i);velocity_x_p(i);velocity_p_kk1(i);position_y_p(i);velocity_y_p(i);velocity_p_kk1(i)];
            %Pre_T_particle(1:6,t,i)=[position_x_p(i);velocity_x_p(i);velocity_p_kk1(i);position_y_p(i);velocity_y_p(i);velocity_p_kk1(i)];
        end
        Pre_T_particle(7,t,:)=1;
        %particle_likehood_after_bias=ones(1,Np);
    else
        %% --------------- evolution of the pre-tracks ----------------
        %% -----------independent partition particle filter------------
        Pre_track_Z=zeros(1,Np);
        Partition_likehood=zeros(1,Np);
        %particle_likehood_after_bias=zeros(1,Np);
        
        for j=1:Np % 6%
            %% === Rao-blackwellisation
            
            %%
              %%%%����ǰ
                                figure(7)
                                colorParticle={'b.','y.','g.','k.';'g^','k^','b^','y^';'bo','ro','mo','go'};
        %                     plot(squeeze(Pre_T_particle(1,t-1,i,:)),squeeze(Pre_T_particle(4,t-1,i,:)),colorParticle{1,i})
                                grid on
                                hold on
                                plot(squeeze(Pre_T_particle(1,t,:)),squeeze(Pre_T_particle(4,t,:)),colorParticle{1,1})
                                title('����ǰ��������') 
            %%
            Pre_T_particle(1:6,t,j)= sample_RB( Pre_T_particle(1:6,t-1,j),T_step,Q_l,Q_n,Q_ln,A_1_t,Q_1_l,q1 ); %����������
            %% 
%                             plot(squeeze(Pre_T_particle(1,t-1,i,:)),squeeze(Pre_T_particle(4,t-1,i,:)),colorParticle{2,i})
                                plot(squeeze(Pre_T_particle(1,t,:)),squeeze(Pre_T_particle(4,t,:)),colorParticle{3,2})
                                xlabel('x��');ylabel('y��');
%                            axis([0,numX,0,numY])
            %%
            Pre_T_life_index(j)=Pre_T_life_index(j)+1;
            Z_x_index=ceil(Pre_T_particle(1,t,j));
            Z_y_index=ceil(Pre_T_particle(4,t,j));
            if Z_x_index<=numX && Z_x_index>0 && Z_y_index<=numY && Z_y_index>0
                Pre_track_Z(j)=Detection_frame(Z_y_index,Z_x_index);
                Pre_T_life_quality(j)=Pre_T_life_quality(j)+Detection_frame(Z_y_index,Z_x_index);
                Pre_T_particle(7,t,j)=Detection_frame(Z_y_index,Z_x_index);
                %% Gaussian likelihood ratio
                Partition_likehood(j)=exp(0.5*(2*Detection_frame(Z_y_index,Z_x_index)*A-A^2));
                %% ��˹����ɢ
%                 for m = 1:numY
%                     for n=1:numX
%                         clean_Frame(m,n) =A*exp(-(((n-Z_x_index)^2/(2)+(m-Z_y_index)^2/(2))*L));
%                     end
%                 end
%                 [index_Gate] = find(10*log(clean_Frame.^2/Sigma_noise)>3);
%                 Partition_likehood(j) = prod(exp(0.5*(2*Detection_frame(index_Gate)*A-A^2)));
                %% Rayleigh likelihood ratio or just likelihood
                %Partition_likehood(j)=raylpdf(Detection_frame(Z_y_index,Z_x_index),sqrt(Sigma_noise+A^2))./raylpdf(Detection_frame(Z_y_index,Z_x_index),Sigma_noise);
%                 Partition_likehood(j)=(1-A/Detection_frame(Z_y_index,Z_x_index))*exp(0.5*(2*Detection_frame(Z_y_index,Z_x_index)*A-A^2)/Sigma_noise^2);
               %% ����ɢ��˹
                   %==================��Ȼ�ȼ���=====================
%                sigma2=Sigma_noise.*ones(1,Np); 
%                mu_h1(j) = Pre_T_particle(7,t,j).*exp(-((Pre_T_particle(1,t,j)-Z_x_index.*Re_x).^2/Re_x+(Pre_T_particle(4,t,j)-Z_y_index.*Re_y).^2/Re_y)*L)+sigma2(j);
%                mu_h0=sigma2;
%  Partition_likehood(j)=((1./mu_h1(j)).*exp(-(1./mu_h1(j)).*Detection_frame(Z_y_index,Z_x_index)))./((1./mu_h0(j)).*exp(-(1./mu_h0(j)).*Detection_frame(Z_y_index,Z_x_index)));
%     
            else
                Partition_likehood(j)=0;
            end
        end
        Partition_likehood=Partition_likehood./sum(Partition_likehood);
        %% === sample index funciton
                                %%%%%%������%%%%%%
%                                 for t=1:5:25
%                                 %%%%����ǰ
%                                 figure(8)
%                                 colorParticle={'b.','y.','g.','k.';'g^','k^','b^','y^';'bo','ro','mo','go'};
%         %                     plot(squeeze(Pre_T_particle(1,t-1,i,:)),squeeze(Pre_T_particle(4,t-1,i,:)),colorParticle{1,i})
%                                 grid on
%                                 hold on
%                                 plot(squeeze(Pre_T_particle(1,t,:)),squeeze(Pre_T_particle(4,t,:)),colorParticle{1,1})
%                                 title('�ز���ǰ��������')
%                                 end
        [index_sample]=Sample_index(Partition_likehood);
        Pre_T_particle(:,t,:)=Pre_T_particle(:,t,index_sample);
        Pre_T_life_quality=Pre_T_life_quality(index_sample);
                                %%%%%%������
                                
%                                 for t=1:5:25
%         %                     plot(squeeze(Pre_T_particle(1,t-1,i,:)),squeeze(Pre_T_particle(4,t-1,i,:)),colorParticle{2,i})
%                                 plot(squeeze(Pre_T_particle(1,t,:)),squeeze(Pre_T_particle(4,t,:)),colorParticle{3,2})
%                                 xlabel('x��');ylabel('y��');
%                                 axis([0,numX,0,numY])
%         %                     title(['��',num2str(i),'��Ŀ���',num2str(t-1),'��',num2str(t),'֡������'])
%         %                     title('����ǰ��������')
%                                 %plot(x_c(3*i-2,t),x_c(3*i-1,t),colorParticle{3,1},'linewidth',2)%,'Ŀ����ʵλ��'
%                                 %legend('����ǰ','������')
%                                 end
        %% === retain the bias of sample: the likehood
        %particle_likehood_after_bias=Partition_likehood(index_sample);
        
    end

end
%% record the estimates
E_target_state=mean(Pre_T_particle(:,:,:),3); %��A�ĵ�dimά�����ֵ
E_target_state_MC([1,2],:)=E_target_state([1,2],:)*Re_x;
E_target_state_MC([4,5],:)=E_target_state([4,5],:)*Re_y;
E_target_state_MC([3,6,7],:)=E_target_state([3,6,7],:);
% %%���θ��ٺ��� E_target_state_MC([1,4],:)=[Re_x,0;0,Re_y]*E_target_state([1,4],:);
%             figure(1)
%             colorP={'b.','r.','g.','k.','','','','','',''};
%             for n=1:Target_number
%                 plot(E_target_state(1,:,n),E_target_state(4,:,n),colorP{n});
%                 if n==1
%                     hold on;xlabel('x��');ylabel('y��');title(['Ŀ��켣��Ŀ�����',num2str(Target_number),'��']);grid on;
%                     %         axis([0,numY,0,numX])
%                 end
%             end






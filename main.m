clc
clear all

M=8; % number of BS antennas
N=16; % number of IRS elements
K=4; % number of users
corr=0; % correlation coefficient between IRS-user channels
BS_pow=43; % BS transmit power in dBm
P_BS=10.^(BS_pow./10).*1e-3; % BS transmit power in Watt
noise_pow=-109; % noise power in dBm/Hz
P1=P_BS; % BS transmit power in Step 1
P2=P_BS; % BS transmit power in Step 2

pilot_time=30:10:80; % total pilot transmitting overhead in time instants
num_case=length(pilot_time);
fb_time=90*ones(1,num_case); % total feedback overhead in time instants
tau1=N*ones(1,num_case); % pilot transmitting overhead in Step 1

mu=4; % QAM modulation rate in bits/sample
bit_allo_mode=1; %bit allocation selection for 1:optimal 2:equal 3:random

bit_ele_ben=fb_time.*mu./(M.*N); % quantization bits for each elements under the benchmark scheme 1
bit_ele_pro=fb_time.*mu./pilot_time; % quantization bits for each elements under the proposed scheme
NMSE_all=zeros(6,num_case);
NMSE_Gk=zeros(3,num_case); % NMSE of cascaded channels under three schemes: 1) proposed scheme; 2) benchmark scheme 1; 3) benchmark scheme 2

num_train=2e4;
num_test=5e3;



for id=1:num_case   

    [G1_all_Tr, coeff_out_Tr]=channel_parameters_downlink(M,N,K,corr,num_train);
    [G1_all, coeff_out]=channel_parameters_downlink(M,N,K,corr,num_test);

    NMSE_all(:,id)=proposed_fb_then_es(M,N,K,noise_pow,P1,P2,G1_all_Tr,coeff_out_Tr,G1_all,coeff_out,bit_ele_pro(:,id),pilot_time(id),tau1(id),mu,bit_allo_mode);
    NMSE_Gk(1,id)=NMSE_all(6,id);
    fprintf('id = %d, NMSE_Gk of proposed with optimal bit allo = %.4f \n',id, NMSE_Gk(1,id));

    % NMSE_Gk(2,id)=proposed_fb_then_es(M,N,K,noise_pow,P1(id),P2(id),G1_all_Tr,coeff_out_Tr,G1_all,coeff_out,bit_ele_pro(id),pilot_time(id),tau1,mu,2);
    % fprintf('id = %d, NMSE_Gk of proposed with equal bit allo = %.4f \n',id, NMSE_Gk(2,id));
    % 
    % NMSE_Gk(3,id)=proposed_fb_then_es(M,N,K,noise_pow,P1(id),P2(id),G1_all_Tr,coeff_out_Tr,G1_all,coeff_out,bit_ele_pro(id),pilot_time(id),tau1,mu,3);
    % fprintf('id = %d, NMSE_Gk of proposed with random bit allo = %.4f \n',id, NMSE_Gk(3,id));

    NMSE_Gk(2,id)=benchmark_es_then_fb(M,N(id),K,noise_pow,BS_pow,G1_all_Tr,coeff_out_Tr,G1_all,coeff_out,bit_ele_ben(id),pilot_time(id),mu);
    fprintf('id = %d, NMSE_Gk of benchmark = %.4f \n',id, NMSE_Gk(2,id));

    NMSE_Gk(3,id)=benchmark_fb_sum(M,N(id),K,noise_pow,BS_pow,G1_all_Tr,coeff_out_Tr,G1_all,coeff_out,fb_time(id),pilot_time(id),mu);
    fprintf('id = %d, NMSE_Gk of benchmark2 = %.4f \n', id, NMSE_Gk(3,id));

    % save('NMSE_Gk.mat','NMSE_Gk');
end
T_tol=fb_time+pilot_time;
figure(1)
for i=1:3
    plot(T_tol,NMSE_Gk(i,:),'o-','LineWidth',1.5);
    grid on,hold on;
    % plot(fb_time,NMSE_all(5,:),'>-','LineWidth',1.5);
    % plot(N,NMSE_Gk(i,:),'>-','LineWidth',1.5);
end
xlabel("Overall Overhead for Both Pilot and Feedback Transmission",'Interpreter', 'latex')
% xlabel('Number of IRS Sub-Surfaces','Interpreter','latex');
ylabel('NMSE of Cascaded Channels','Interpreter', 'latex');
% labels=["Proposed Scheme","Benchmark Scheme 1","Benchmark Scheme 2"];
% legend(labels,'Interpreter', 'latex');

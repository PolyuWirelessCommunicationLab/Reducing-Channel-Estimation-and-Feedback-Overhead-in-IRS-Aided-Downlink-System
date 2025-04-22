clear
M=12;
N=16;
K=11;
corr=0.01;
BS_pow=43;%dBm
P_BS=10.^(BS_pow./10).*1e-3; %Watt
noise_pow=-109;
tau2_min=ceil(N.*(M-1)./K);


fb_time=78:12:132;
num_case=length(fb_time);
pilot_time=48*ones(1,num_case);

tau1=N*ones(1,num_case);
tau2=pilot_time-tau1;
% tau2=pilot_time./(1+tau_ratio);
% tau1=pilot_time-tau2;

P_sum=P_BS;
P1=P_BS;
P2=P_BS;

mu=4;
pilot_irs_mode=2;
bit_allo_mode=1; %1:optimal 2:equal 3:random

bit_ele_ben=fb_time.*mu./(M*N);
bit_ele_ben2=bit_ele_ben;
bit_ele_pro=fb_time.*mu./pilot_time;

NMSE_all=zeros(6,num_case);
NMSE_Gk=zeros(4,num_case);

num_train=1e4;
num_test=5e3;

[G1_all_Tr, coeff_out_Tr]=channel_parameters_downlink(M,N,K,corr,num_train);
[G1_all, coeff_out]=channel_parameters_downlink(M,N,K,corr,num_test);

for id=1:num_case        

    NMSE_Gk(1,id)=proposed_fb_then_es(M,N,K,noise_pow,P1(id),P2(id),G1_all_Tr,coeff_out_Tr,G1_all,coeff_out,bit_ele_pro(:,id),pilot_time(id),tau1(id),mu,pilot_irs_mode,bit_allo_mode);
    fprintf('id = %d, NMSE_Gk of proposed with optimal bit allo = %.4f \n',id, NMSE_Gk(1,id));

    % NMSE_Gk(2,id)=proposed_fb_then_es(M,N,K,noise_pow,P1(id),P2(id),G1_all_Tr,coeff_out_Tr,G1_all,coeff_out,bit_ele_pro(id),pilot_time(id),tau1,mu,pilot_irs_mode,2);
    % fprintf('id = %d, NMSE_Gk of proposed with equal bit allo = %.4f \n',id, NMSE_Gk(2,id));

    % NMSE_Gk(3,id)=proposed_fb_then_es(M,N,K,noise_pow,P1(id),P2(id),G1_all_Tr,coeff_out_Tr,G1_all,coeff_out,bit_ele_pro(id),pilot_time(id),tau1,mu,pilot_irs_mode,3);
    % fprintf('id = %d, NMSE_Gk of proposed with random bit allo = %.4f \n',id, NMSE_Gk(3,id));

    % NMSE_Gk(3,id)=benchmark_es_then_fb(M,N,K,noise_pow,BS_pow./pilot_time(id),G1_all_Tr,coeff_out_Tr,G1_all,coeff_out,bit_ele_ben(id),pilot_time(id),mu);
    % fprintf('id = %d, NMSE_Gk of benchmark = %.4f \n',id, NMSE_Gk(3,id));

    % NMSE_Gk(4,id)=benchmark_fb_sum(M,N,K,noise_pow,BS_pow./pilot_time(id),G1_all_Tr,coeff_out_Tr,G1_all,coeff_out,fb_time(id),pilot_time(id),mu);
    % fprintf('id = %d, NMSE_Gk of benchmark2 = %.4f \n', id, NMSE_Gk(4,id));

    % save('NMSE_Gk.mat','NMSE_Gk');
end
T_tol=fb_time+pilot_time;
figure(1)
for i=[1]
    % plot(fb_time,NMSE_all(4,:),'o--','LineWidth',1.5);
    grid on,hold on;
    % plot(fb_time,NMSE_all(5,:),'>-','LineWidth',1.5);
    plot(pow_ratio,NMSE_Gk(i,:),'o-','LineWidth',1.5);
end
% xlabel('$$\bar{p}:\tilde{p}$$','Interpreter', 'latex'); 
% xlabel("Overall Overhead for Both Pilot and Feedback Transmission")
ylabel('NMSE of Cascaded Channels','Interpreter', 'latex');
% xlim([0.1 30])
% labels=["Proposed Scheme","Benchmark Scheme 1","Benchmark Scheme 2"];
labels=["$$\tau_1:\tau_2=1/5$$","$$\tau_1:\tau_2=1$$","$$\tau_1:\tau_2=3$$"];
legend(labels,'Interpreter', 'latex');
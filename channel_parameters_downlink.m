function [G1_all, coeff_out]=channel_parameters_downlink(M,N,K,corr,N_para)
%% ==========================
% generating channels
% G1_all : user 1's overall reflecting channel
% coeff_out : correlation coefficients of the channels of other users to user 1
%% ===========================
c_irs=corr; % BS-IRS correlation
c_bs=0.01; % BS correlation
d0=0.01;
d_bs=d0*100^(-2.2); % pathloss between IRS and BS
alph_bi=10;

C_irs=C_gen(N,c_irs); % generating of square root correlation matrix
C_bs=C_gen(M,c_bs); % generating of square root correlation matrix
c_irs_ur=corr*ones(1,K); % IRS-user correlation
C_irs_ur=zeros(K*N,N);
for k=1:K
    C_irs_ur((k-1)*N+1:k*N,:)=C_gen(N,c_irs_ur(k));
end

G1_all=zeros(M*N_para,N);
coeff=zeros((K-1)*N,N_para);

for i=1:N_para   
    rng(i,"twister");  
    dist_IRS_UE=generate_user_positions(K);
    d_ue=d0*dist_IRS_UE.^(-2.1); % pathloss between IRS and users
    % reflecting chanel of user 1
    s_upa=UPA_gen(N);
    s_ula=ULA_gen(M); % column vector
    R_det=s_ula*s_upa'; % the los component between BS and IRS
    H_bs_hat=sqrt(d_bs/2)*(randn(M,N)+1j*randn(M,N));
    
    % Racian
    H_bs=sqrt(d_bs)*sqrt(alph_bi/(1+alph_bi))*R_det+sqrt(1/(1+alph_bi))*C_bs*H_bs_hat*C_irs;
    % Rayleigh
    % H_bs=C_bs*H_bs_hat*C_irs;

    C_irs_1=C_irs_ur(1:N,:);
    H_ue1_hat=sqrt(d_ue(1)/2)*(randn(N,1)+1j*randn(N,1));
    H_ue1=C_irs_1*H_ue1_hat;
    H1_mtx=H_bs*diag(H_ue1);
    G1_all((i-1)*M+1:i*M,:)=H1_mtx;
    % reflecting channel of users 2-K; coefficients
      for k=2:K
         % generating matrix from user 2 to K.
         C_irs_k=C_irs_ur((k-1)*N+1:k*N,:);
         H_uek_hat=sqrt(d_ue(k)/2)*(randn(N,1)+1j*randn(N,1));
         H_uek=C_irs_k*H_uek_hat;
         coeff((k-2)*N+1:(k-1)*N,i)=H_uek./H_ue1;
      end
end
coeff_out=coeff;
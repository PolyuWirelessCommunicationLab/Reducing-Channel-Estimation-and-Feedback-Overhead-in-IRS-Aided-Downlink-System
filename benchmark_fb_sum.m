function NMSE_Gk=benchmark_fb_sum(M,N,K,noise_pow,BS_pow,G1_all_Tr,coeff_out_Tr,G1_all,coeff_out,fb_time,pilot_length,mu)

tau=pilot_length; %time slots
P1=10^(BS_pow/10)*1e-3;  %transmit power 
Phi1=dftmtx(N);
Phi2=dftmtx(max(N,pilot_length-N));
x1=zeros(tau,M);
x1(1:N,:)=ones(N,M);
x1(N+1:end,:)=1/sqrt(2)*(randn(pilot_length-N,M)+1j*randn(pilot_length-N,M));
A1=zeros(tau,M*N);
for i=1:tau
    if i<=N
        Ai=Phi1(i,1:N)*(sqrt(P1)*kron(eye(N),x1(i,:)));
    else
        Ai=Phi2(i-N,1:N)*(sqrt(P1)*kron(eye(N),x1(i,:)));
    end
    A1(i,:)=Ai;
end
num_train=size(coeff_out_Tr,2);
num_test=size(coeff_out,2);

Gk_set_Tr=zeros(M*N*K,num_train);
Y_pure_Tr=zeros(tau*K,num_train);
Alpha_set_Tr=zeros(K*N,num_train);
Cov_Gk=zeros(M*N*K,M*N);
Cov_Alpha=zeros(K*N,N);
for ii=1:num_train
    coeff_ii=reshape(coeff_out_Tr(:,ii),[N,K-1]);
    LAMDA=[ones(1,N);(coeff_ii).']; % reflecting channel ratio
    G1=G1_all_Tr((ii-1)*M+1:ii*M,:);    
    Gk=repmat(G1,[K,1]).*repelem(LAMDA,M,1); % userk 's reflecting channel
    for k=1:K
        Gk_k=Gk(M*(k-1)+1:M*k,:);        
        Alpha_ki((k-1)*N+1:k*N)=(sqrt(P1)*x1(1,:)*Gk_k).';        
        Y_pure_Tr((k-1)*tau+1:k*tau,ii)=A1*Gk_k(:);
        Gk_set_Tr(M*N*(k-1)+1:M*N*k,ii)=Gk_k(:);
        Cov_Gk((k-1)*M*N+1:k*M*N,:)=Cov_Gk((k-1)*M*N+1:k*M*N,:)+(Gk_k(:))*(Gk_k(:))';
        Cov_Alpha((k-1)*N+1:k*N,:)=Cov_Alpha((k-1)*N+1:k*N,:)+(Alpha_ki((k-1)*N+1:k*N))'*Alpha_ki((k-1)*N+1:k*N);
    end
    Alpha_set_Tr(:,ii)=Alpha_ki;
end
Cov_Gk=Cov_Gk./num_train;
Cov_Alpha=Cov_Alpha./num_train;


sigma_z=sqrt(10^(noise_pow/10)*1e-3); 
noise=zeros(tau*K,num_train);
for k=1:K
    noise((k-1)*tau+1:k*tau,:)=1/sqrt(2)*(randn(tau,num_train)+1j*randn(tau,num_train));
end
Y_set_Tr=Y_pure_Tr+sigma_z.*noise;  

Gk_set=zeros(M*N*K,num_test);
Y_pure=zeros(tau*K,num_test);
gk_set=zeros(K*N,num_test);
Alpha_set=zeros(K*N,num_test);

for ii=1:num_test
    coeff_ii=reshape(coeff_out(:,ii),[N,K-1]);
    LAMDA=[ones(1,N);(coeff_ii).']; % reflecting channel ratio
    G1=G1_all((ii-1)*M+1:ii*M,:);      
    Gk=repmat(G1,[K,1]).*repelem(LAMDA,M,1); % userk 's reflecting channel
    for k=1:K
        Gk_k=Gk(M*(k-1)+1:M*k,:);        
        Alpha_ki((k-1)*N+1:k*N)=(sqrt(P1)*ones(1,M)*Gk_k).';        
        Y_pure((k-1)*tau+1:k*tau,ii)=A1*Gk_k(:);
        Gk_set(M*N*(k-1)+1:M*N*k,ii)=Gk_k(:);
        % gk_k_d=sum(Gk_k,1);
        % gk_set((k-1)*N+1:k*N,ii)=gk_k_d.';
    end
    Alpha_set(:,ii)=Alpha_ki;
end
noise=zeros(tau*K,num_test);
for k=1:K
    noise((k-1)*tau+1:k*tau,:)=1/sqrt(2)*(randn(tau,num_test)+1j*randn(tau,num_test));
end
Y_set=Y_pure+sigma_z.*noise;  
%% Estimation and Quantization
Gk_es_set=zeros(M*N*K,num_test);
Gk_es_set_Tr=zeros(M*N*K,num_train);
gk_es_set=zeros(N*K,num_test);
Alpha_es_set_Tr=zeros(K*N,num_train);
Alpha_es_set=zeros(K*N,num_test);
for k=1:K
    Cov_Gk_k=Cov_Gk((k-1)*M*N+1:k*M*N,:);      
    Y1k_Tr_2=Y_set_Tr((k-1)*tau+N+1:k*tau,:);
    Y1k_2=Y_set((k-1)*tau+N+1:k*tau,:);
    Gk_es_Tr=Cov_Gk_k*A1'/(A1*Cov_Gk_k*A1'+sigma_z^2*eye(pilot_length))*Y_set_Tr((k-1)*tau+1:k*tau,:);
    Gk_es_k=Cov_Gk_k*A1'/(A1*Cov_Gk_k*A1'+sigma_z^2*eye(pilot_length))*Y_set((k-1)*tau+1:k*tau,:);
    Gk_es_set((k-1)*M*N+1:k*M*N,:)=Gk_es_k;
    Gk_es_set_Tr((k-1)*M*N+1:k*M*N,:)=Gk_es_Tr;

    Cov_Alphak=Cov_Alpha((k-1)*N+1:k*N,:);
    Y1k_Tr_1=Y_set_Tr((k-1)*tau+1:(k-1)*tau+N,:);
    Alphak_es_Tr=Cov_Alphak*Phi1'/(Phi1*Cov_Alphak*Phi1'+sigma_z^2*eye(N))*Y1k_Tr_1;
    Alpha_es_set_Tr((k-1)*N+1:k*N,:)=Alphak_es_Tr;
    Y1k_1=Y_set((k-1)*tau+1:(k-1)*tau+N,:);
    Alphak_es=Cov_Alphak*Phi1'/(Phi1*Cov_Alphak*Phi1'+sigma_z^2*eye(N))*Y1k_1;
    Alpha_es_set((k-1)*N+1:k*N,:)=Alphak_es;
   
    % for d=1:N
    %     Gk_es_d=Gk_es_k((d-1)*M+1:d*M,:);
    %     gk_es_d=sum(Gk_es_d,1);
    %     gk_es_set((k-1)*N+d,:)=gk_es_d;
    % end
end
MSE_Gk_es=sum(abs(Gk_es_set-Gk_set).^2,1);
pow_Gk=sum(abs(Gk_set).^2,1);
NMSE_Gk_es=mean(MSE_Gk_es)./mean(pow_Gk);

MSE_Alpha_es=sum(abs(Alpha_es_set-Alpha_set).^2,1);
pow_Alpha=sum(abs(Alpha_set).^2,1);
NMSE_Alpha_es=mean(MSE_Alpha_es)./mean(pow_Alpha);
NMSE_all(1)=NMSE_Alpha_es;

%quantization of all users' alphas

bit_ele=fb_time/(M*N/K+N)*mu;

fb_ratio=0.45;
bit_sum=fb_time*fb_ratio/N*mu;
bits_sum=nearest_integers_avg(bit_ele,N);
Alpha_quan_set=zeros(K*N,num_test);
for k=1:K
    Alphak_Train=Alpha_es_set_Tr((k-1)*N+1:k*N,:); 
    Alphak_es=Alpha_es_set((k-1)*N+1:k*N,:);
    Alpha_quan=zeros(N,num_test);
    parfor d=1:N
        [codebook_Alpha,~]=VQ_Lloyd(Alphak_Train(d,:),bits_sum(d));
        Alpha_quan_d=Quantizer(Alphak_es(d,:),codebook_Alpha);      
        Alpha_quan(d,:)=Alpha_quan_d;
    end
    Alpha_quan_set((k-1)*N+1:k*N,:)=Alpha_quan;  
end



% for k=1:K
%     Gk_k_Train=Gk_es_set_Tr((k-1)*M*N+1:k*M*N,:); 
%     Gk_es_k=Gk_es_set((k-1)*M*N+1:k*M*N,:);
%     gk_quan=zeros(N,num_test);
%     Gk_Train=zeros(N,M,num_train);
%     Gk_es=zeros(N,M,num_test);
%     for d=1:N
%         Gk_Train(d,:,:)=Gk_k_Train((d-1)*M+1:d*M,:);
%         Gk_es(d,:,:)=Gk_es_k((d-1)*M+1:d*M,:);
%     end
%     parfor d=1:N
%         Gk_d_Train=squeeze(Gk_Train(d,:,:)); %M*num_train
%         Gk_es_d=squeeze(Gk_es(d,:,:));
%         gk_d_Train=sum(Gk_d_Train,1); %sum of g_{k,d}, 1*num_train           
%         gk_es_d=sum(Gk_es_d,1); %sum of g_es_{k,d}, 1*num_test
% 
%         [codebook_gkd,~]=VQ_Lloyd(gk_d_Train,bits_sum(d));
%         gk_quan_d=Quantizer(gk_es_d,codebook_gkd);      
%         gk_quan(d,:)=gk_quan_d;
%     end
%     gk_quan_set((k-1)*N+1:k*N,:)=gk_quan;  
% end
MSE_Alpha=sum(abs(Alpha_quan_set-Alpha_set).^2,1);
pow_Alpha=sum(abs(Alpha_set).^2,1);
NMSE_Alpha_quan=mean(MSE_Alpha)./mean(pow_Alpha);
NMSE_all(2)=NMSE_Alpha_quan;

%index of reference users
ref_idx=zeros(N,num_test);
sk=zeros(K,num_test);
for ii=1:num_test
    Alpha_quan=reshape(Alpha_quan_set(:,ii),N,K);       
    for d=1:N
        [~,ref_idx(d,ii)]=max(abs(Alpha_quan(d,:)));
        % ref_idx(d,ii)=1;
    end
    for k=1:K
        sk(k,ii)=sum(ref_idx(:,ii)==k);        
    end
end    


%quantization of reference users's channels
Gref_set=zeros(M*N,num_test);
Gref_quan=zeros(M*N,num_test); 

bit_ref=fb_time*(1-fb_ratio)./(M*N/K)*mu;
T_fb=(bit_sum*N+bit_ref*M*N/K)/mu;

bits_ref=nearest_integers_avg(bit_ele,M*N);
for k=1:K        
    Gk_k_Train=Gk_es_set_Tr((k-1)*M*N+1:k*M*N,:); 
    Gk_es_k=Gk_es_set((k-1)*M*N+1:k*M*N,:);
    Gk_k_real=Gk_set((k-1)*M*N+1:k*M*N,:);

    Gk_k_quan=zeros(M*N,num_test);
    parfor n=1:N*M
        Gkn_Train=Gk_k_Train(n,:);
        [codebook_Gk,~]=VQ_Lloyd(Gkn_Train,bits_ref(n));
        Gkn_seg=Gk_es_k(n,:);
        Gkn_quan=Quantizer(Gkn_seg,codebook_Gk);
        Gk_k_quan(n,:)=Gkn_quan;
    end
    
    for d=1:N
        % Gk_k_d_Tr=Gk_k_Train((d-1)*M+1:d*M,:);
        % Gk_es_k_d=Gk_es_k((d-1)*M+1:d*M,:);
        % parfor m=1:M
        %     [codebook_Gk,~]=VQ_Lloyd(Gk_k_d_Tr(m,:),bit_allo(m));
        %     Gk_k_d_quan(m,:)=Quantizer(Gk_es_k_d(m,:),codebook_Gk); 
        % end
        for ii=1:num_test
            if k==ref_idx(d,ii)
                Gref_quan((d-1)*M+1:d*M,ii)=Gk_k_quan((d-1)*M+1:d*M,ii);
                Gref_set((d-1)*M+1:d*M,ii)=Gk_k_real((d-1)*M+1:d*M,ii);                
            end
        end
    end    
end
g_ref_quan=zeros(N,num_test);
for d=1:N
    g_ref_quan(d,:)=sum(Gref_quan((d-1)*M+1:d*M,:));
end
MSE_Gref=sum(abs(Gref_quan-Gref_set).^2,2);
pow_Gref=sum(abs(Gref_set).^2,2);
NMSE_Gref=mean(MSE_Gref)./mean(pow_Gref);
NMSE_all(3)=NMSE_Gref;

%calculation of coefficients
coeff_quan_set=zeros(K*N,num_test);
coeff_new_set=zeros(K*N,num_test);
parfor ii=1:num_test
    coeff_quan=zeros(N,K);
    coeff_new=zeros(N,K);
    Alpha_quan=reshape(Alpha_quan_set(:,ii),N,K);
    Alpha=reshape(Alpha_set(:,ii),N,K);
    for d=1:N           
        coeff_quan(d,:)=Alpha_quan(d,:)./Alpha_quan(d,ref_idx(d,ii));
        coeff_new(d,:)=Alpha(d,:)./Alpha(d,ref_idx(d,ii));
    end
    coeff_quan_set(:,ii)=coeff_quan(:);
    coeff_new_set(:,ii)=coeff_new(:);
end
%recover of all gk's
Gk_quan_set=zeros(K*M*N,num_test);
% Gk_new_set=zeros(K*M*N,num_test);
for k=1:K        
    Gk_quan_set(M*N*(k-1)+1:M*N*k,:)=Gref_quan.*repelem(coeff_quan_set((k-1)*N+1:k*N,:),M,1);
    % Gk_new_set(M*N*(k-1)+1:M*N*k,:)=Gref_set.*repelem(coeff_new_set((k-1)*N+1:k*N,:),M,1);        
end
MSE_coeff=sum(abs(coeff_quan_set-coeff_new_set).^2,1);
pow_coeff=sum(abs(coeff_new_set).^2,1);
NMSE_coeff=mean(MSE_coeff)./mean(pow_coeff);
NMSE_all(4)=NMSE_coeff;

MSE_Gk=sum(abs(Gk_quan_set-Gk_set).^2,1);
pow_Gk=sum(abs(Gk_set).^2,1);
NMSE_Gk=mean(MSE_Gk)./mean(pow_Gk);   
NMSE_all(5)=NMSE_Gk;

T_tol=pilot_length+fb_time;

end 

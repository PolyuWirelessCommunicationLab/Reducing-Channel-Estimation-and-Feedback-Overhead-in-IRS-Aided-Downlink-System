function NMSE_Gk=benchmark_es_then_fb(M,N,K,noise_pow,BS_pow,G1_all_Tr,coeff_out_Tr,G1_all,coeff_out,bit_ele,pilot_length,mu)
%%
% Benchmark Scheme 1
%%
tau=pilot_length; %time slots
P1=10^(BS_pow/10)*1e-3;  %transmit power 
Phi1=dftmtx(tau)/sqrt(tau);

x1=1/sqrt(2)*(randn(tau,M)+1j*randn(tau,M));
A1=zeros(tau,M*N);
for i=1:tau
    Ai=Phi1(i,1:N)*(sqrt(P1)*kron(eye(N),x1(i,:)));
    A1(i,:)=Ai;
end
num_train=size(coeff_out_Tr,2);
num_test=size(coeff_out,2);

Gk_set_Tr=zeros(M*N*K,num_train);
Cov_Gk=zeros(M*N*K,M*N);
mean_Gk=zeros(K*M*N,1);
Y_pure_Tr=zeros(tau*K,num_train);
for ii=1:num_train
    coeff_ii=reshape(coeff_out_Tr(:,ii),[N,K-1]);
    LAMDA=[ones(1,N);(coeff_ii).']; % reflecting channel ratio
    G1=G1_all_Tr((ii-1)*M+1:ii*M,:);    
    Gk=repmat(G1,[K,1]).*repelem(LAMDA,M,1); % userk 's reflecting channel
    for k=1:K
        Gk_k=Gk(M*(k-1)+1:M*k,:);        
        Gk_set_Tr(M*N*(k-1)+1:M*N*k,ii)=Gk_k(:);
        Y_pure_Tr((k-1)*tau+1:k*tau,ii)=A1*Gk_k(:);
        Cov_Gk((k-1)*M*N+1:k*M*N,:)=Cov_Gk((k-1)*M*N+1:k*M*N,:)+(Gk_k(:))*(Gk_k(:))';
    end
end
Cov_Gk=Cov_Gk./num_train;
sigma_z=sqrt(10^(noise_pow/10)*1e-3); 
noise=zeros(tau*K,num_train);
for k=1:K
    noise((k-1)*tau+1:k*tau,:)=1/sqrt(2)*(randn(tau,num_train)+1j*randn(tau,num_train));
end
Y_set_Tr=Y_pure_Tr+sigma_z.*noise;

Gk_es_set_Tr=zeros(M*N*K,num_train);
for k=1:K
    Cov_Gk_k=Cov_Gk((k-1)*M*N+1:k*M*N,:);    
    Gk_es_k_Tr=Cov_Gk_k*A1'/(A1*Cov_Gk_k*A1'+sigma_z^2*eye(tau))*Y_set_Tr((k-1)*tau+1:k*tau,:);       
    Gk_es_set_Tr((k-1)*M*N+1:k*M*N,:)=Gk_es_k_Tr;
    Gk_k_Tr=Gk_set_Tr((k-1)*M*N+1:k*M*N,:);
    MSE_Gk_Tr=sum(abs(Gk_k_Tr-Gk_es_k_Tr).^2,1);
    pow_Gk_Tr=sum(abs(Gk_k_Tr).^2,1);
    NMSE_GK_es_Tr=mean(MSE_Gk_Tr)./mean(pow_Gk_Tr);
end


Gk_set=zeros(M*N*K,num_test);
Y_pure=zeros(tau*K,num_test);
Cov_Gk_test=zeros(K*M*N,M*N);
for ii=1:num_test
    coeff_ii=reshape(coeff_out(:,ii),[N,K-1]);
    LAMDA=[ones(1,N);(coeff_ii).']; % reflecting channel ratio
    G1=G1_all((ii-1)*M+1:ii*M,:);      
    Gk=repmat(G1,[K,1]).*repelem(LAMDA,M,1); % userk 's reflecting channel
    for k=1:K
        Gk_k=Gk(M*(k-1)+1:M*k,:);        
        Y_pure((k-1)*tau+1:k*tau,ii)=A1*Gk_k(:);
        Gk_set(M*N*(k-1)+1:M*N*k,ii)=Gk_k(:);
        Cov_Gk_test((k-1)*M*N+1:k*M*N,:)=Cov_Gk_test((k-1)*M*N+1:k*M*N,:)+(Gk_k(:))*(Gk_k(:))';
    end
end
Cov_Gk_test=Cov_Gk_test./num_test;
sigma_z=sqrt(10^(noise_pow/10)*1e-3); 
noise=zeros(tau*K,num_test);
for k=1:K
    noise((k-1)*tau+1:k*tau,:)=1/sqrt(2)*(randn(tau,num_test)+1j*randn(tau,num_test));
end
Y_set=Y_pure+sigma_z.*noise;

bit_allo=nearest_integers_avg(bit_ele,M*N);
Gk_quan_set=zeros(K*M*N,num_test);
for k=1:K
    Gk_k_Train=Gk_es_set_Tr(M*N*(k-1)+1:M*N*k,:); 
    Cov_Gk_k=Cov_Gk((k-1)*M*N+1:k*M*N,:);

    Gk_quan=zeros(M*N,num_test);  
    Gk_es_k=Cov_Gk_k*A1'/(A1*Cov_Gk_k*A1'+sigma_z^2*eye(tau))*Y_set((k-1)*tau+1:k*tau,:);
    Gk_k=Gk_set((k-1)*M*N+1:k*M*N,:);
    MSE_Gk_es=sum(abs(Gk_k-Gk_es_k).^2,1);
    pow_Gk=sum(abs(Gk_k).^2,1);
    NMSE_Gk_es=mean(MSE_Gk_es)./mean(pow_Gk);

    parfor n=1:M*N
        if bit_allo(n)~=0
            Gkn_Train=Gk_k_Train(n,:);
            [codebook_Gkn,~]=VQ_Lloyd(Gkn_Train,bit_allo(n));
            Gkn_seg=Gk_es_k(n,:);  
            Gkn_quan=Quantizer(Gkn_seg,codebook_Gkn);
            Gk_quan(n,:)=Gkn_quan;
        end
    end 
    Gk_quan_set(M*N*(k-1)+1:M*N*k,:)=Gk_quan;   
end    
T_tol=tau+bit_ele*M*N/mu;
MSE_Gk=sum(abs(Gk_quan_set-Gk_set).^2,1);
pow_Gk=sum(abs(Gk_set).^2,1);
NMSE_Gk=mean(MSE_Gk)./mean(pow_Gk);
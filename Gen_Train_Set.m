function [Y_set_Tr,Cov_Y,Y_set_Test,Alpha_set_test,Y2_pure_set]=Gen_Train_Set(G1_all_Tr,coeff_out_Tr,G1_all,coeff_out,M,N,K,tau1,tau2,tau_tol,P1,P2,noise_pow,Phi1,x1,x2,mode)
%%
% Generate date set of training and testing for Lloyd algorithm 
% Y_set_Tr: training set of received signals
% Cov_Y: covariance matrix of received signals
% Y_set_Test: test set of received signals
% Alpha_set_test: test set of \alpha_{k,d}'s
% Y2_pure_set: training set of received signals without noise in Step 2 
%%
num_train=size(coeff_out_Tr,2);
num_test=size(coeff_out,2);
% Training Set
Alpha_set_Tr=zeros(N*K,num_train);
Cov_Alpha=zeros(N*K,N);
Y1_pure_set=zeros(tau1*K,num_test);
Y2_pure_set=zeros(tau2*K,num_train);

for ii=1:num_train
    coeff_ii=reshape(coeff_out_Tr(:,ii),[N,K-1]);
    LAMDA=[ones(1,N);(coeff_ii).']; % reflecting channel ratio
    G1=G1_all_Tr((ii-1)*M+1:ii*M,:);    
    Gk=repmat(G1,[K,1]).*repelem(LAMDA,M,1); % userk 's reflecting channel
    Alpha=(sqrt(P1)*kron(eye(K),x1.')*Gk).';
    Y1_pure=Phi1*Alpha;
    Y1_pure_set(:,ii)=Y1_pure(:);
    Alpha_set_Tr(:,ii)=Alpha(:);
    for k=1:K
       Cov_Alpha((k-1)*N+1:k*N,:)=Cov_Alpha((k-1)*N+1:k*N,:)+Alpha_set_Tr((k-1)*N+1:k*N,ii)*(Alpha_set_Tr((k-1)*N+1:k*N,ii))';
    end    
    THETA=Gen_THETA(K,M,N,tau2,P2,LAMDA,x2);
    OMEGA=[THETA;kron(eye(N),sqrt(P1)*x1.')]; 
    Y2_pure_set(:,ii)=OMEGA(1:tau2*K,:)*G1(:);
end

%AWGN
sigma_z1=sqrt(10^(noise_pow/10)*1e-3);
noise1=zeros(tau1*K,num_train);
for k=1:K
    noise1((k-1)*tau1+1:k*tau1,:)=1/sqrt(2)*(randn(tau1,num_train)+1j*randn(tau1,num_train));
end
Y1_set_Tr=Y1_pure_set+sigma_z1.*noise1;
Cov_Alpha=Cov_Alpha./num_train;

Alpha_es_set_Tr=zeros(K*N,num_train);
for k=1:K
    Cov_Alphak=Cov_Alpha((k-1)*N+1:k*N,:);
    Y1k=Y1_set_Tr((k-1)*tau1+1:k*tau1,:);
    Alphak_es=Cov_Alphak*Phi1'/(Phi1*Cov_Alphak*Phi1'+sigma_z1^2*eye(tau1))*Y1k;
    Alpha_es_set_Tr((k-1)*N+1:k*N,:)=Alphak_es;   
end


sigma_z2=sqrt(10^(noise_pow/10)*1e-3);
noise2=zeros(tau2*K,num_train);
for k=1:K
    noise2((k-1)*tau2+1:k*tau2,:)=1/sqrt(2)*(randn(tau2,num_train)+1j*randn(tau2,num_train));
end
Y2_set=Y2_pure_set+sigma_z2.*noise2;

Y_set_Tr=zeros(K*tau_tol,num_train);
Cov_Y=zeros(K*tau_tol,tau_tol);
for k=1:K
    for ii=1:num_train
        Y_set_Tr((k-1)*tau_tol+1:k*tau_tol,ii)=[Alpha_es_set_Tr((k-1)*N+1:k*N,ii);Y2_set(k:K:K*(tau2-1)+k,ii)];
        Cov_Y((k-1)*tau_tol+1:k*tau_tol,:)=Cov_Y((k-1)*tau_tol+1:k*tau_tol,:)+Y_set_Tr((k-1)*tau_tol+1:k*tau_tol,ii)*(Y_set_Tr((k-1)*tau_tol+1:k*tau_tol,ii))';
    end
end
Cov_Y=Cov_Y./num_train;

%Test set
Y1_pure_set=zeros(tau1*K,num_test);
Alpha_set_test=zeros(N*K,num_test);
Y2_pure_set=zeros(tau2*K,num_test);

for ii=1:num_test
    coeff_ii=reshape(coeff_out(:,ii),[N,K-1]);
    LAMDA=[ones(1,N);(coeff_ii).']; % reflecting channel ratio
    G1=G1_all((ii-1)*M+1:ii*M,:);    
    Gk=repmat(G1,[K,1]).*repelem(LAMDA,M,1); % userk 's reflecting channel
    Alpha=(sqrt(P1)*kron(eye(K),x1.')*Gk).';
    Y1_pure=Phi1*Alpha;
    Y1_pure_set(:,ii)=Y1_pure(:);
    Alpha_set_test(:,ii)=Alpha(:);
    THETA=Gen_THETA(K,M,N,tau2,P2,LAMDA,mode,x2);
    OMEGA=[THETA;kron(eye(N),sqrt(P1)*x1.')]; 
    Y2_pure_set(:,ii)=OMEGA(1:tau2*K,:)*G1(:);
end

%AWGN
sigma_z1=sqrt(10^(noise_pow/10)*1e-3);
noise1=zeros(tau1*K,num_test);
for k=1:K
    noise1((k-1)*tau1+1:k*tau1,:)=1/sqrt(2)*(randn(tau1,num_test)+1j*randn(tau1,num_test));
end
Y1_set=Y1_pure_set+sigma_z1.*noise1;  

Alpha_es_set=zeros(K*N,num_test);
for k=1:K
    Cov_Alphak=Cov_Alpha((k-1)*N+1:k*N,:);
    Y1k=Y1_set((k-1)*tau1+1:k*tau1,:);
    Alphak_es=Cov_Alphak*Phi1'/(Phi1*Cov_Alphak*Phi1'+sigma_z1^2*eye(tau1))*Y1k;
    Alpha_es_set((k-1)*N+1:k*N,:)=Alphak_es;
end

sigma_z2=sqrt(10^(noise_pow/10)*1e-3);
noise2=zeros(tau2*K,num_test);
for k=1:K
    noise2((k-1)*tau2+1:k*tau2,:)=1/sqrt(2)*(randn(tau2,num_test)+1j*randn(tau2,num_test));
end
Y2_set=Y2_pure_set+sigma_z2.*noise2;

Y_set_Test=zeros(K*tau_tol,num_test);
for k=1:K
    for ii=1:num_test
        Y_set_Test((k-1)*tau_tol+1:k*tau_tol,ii)=[Alpha_es_set((k-1)*N+1:k*N,ii);Y2_set(k:K:K*(tau2-1)+k,ii)];
    end
end


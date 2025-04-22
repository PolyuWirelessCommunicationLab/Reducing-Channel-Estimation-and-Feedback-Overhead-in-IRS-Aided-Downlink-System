function NMSE_Y=proposed_fb_then_es(M,N,K,noise_pow,P1,P2,G1_all_Tr,coeff_out_Tr,G1_all,coeff_out,bit_ele,pilot_length,tau1,mu,pilot_irs_mode,bit_allo_mode)

tau2=pilot_length-tau1;
% tau2=ceil(N*(M-1)/K);
% tau1=pilot_length-tau2; %time slots
x1=ones(M,1); % pilot signal
F=dftmtx(max(tau1,N));
Phi1=F(1:tau1,1:N);
% tau2=M*N-tau1;
tau_tol=N+tau2;

x2=1/sqrt(2)*(randn(tau2,M)+1j*randn(tau2,M));
% x2(1:3,:)=x2(1:3,:)/4;
num_test=size(coeff_out,2);

[Y_set_Tr,Cov_Y,Y_set_Test,Alpha_set,Y2_pure_set]=Gen_Train_Set(G1_all_Tr,coeff_out_Tr,G1_all,coeff_out,M,N,K,tau1,tau2,tau_tol,P1,P2,noise_pow,Phi1,x1,x2,pilot_irs_mode);

%total overhead
Bits_tol=bit_ele.*tau_tol;
T_tol=pilot_length+max(Bits_tol/mu);

bit_allo=zeros(tau_tol,K);
%bit allocation strategy 1: optimized 2: equal 3: random
switch bit_allo_mode
    case 1
        bit_allo=bit_allo_allrev(tau_tol, N, K, Cov_Y, Bits_tol);
    case 2
        bits=zeros(tau_tol,1);
        % bits=(nearest_integers_avg(bit_ele,tau_tol)).';
        bits(1:N)=nearest_integers_avg(bit_ele(1),N);
        bits(N+1:end)=nearest_integers_avg(sum(bit_ele),tau2);
        bit_allo=repmat(bits,1,K);
    case 3
        bit_allo=zeros(tau_tol,K);
        for k=1:K
            bit_allo(:,k)=generate_vector(tau_tol,Bits_tol);
        end
end
codebook_Y=cell(K,tau_tol);
%Codebook
for k=1:K
    Cov_Yk=Cov_Y((k-1)*tau_tol+1:k*tau_tol,:);
    pow_Yk=diag(Cov_Yk);
    scale_fac=mean(pow_Yk(1:N))/mean(pow_Yk(N+1:tau_tol));
    [Pk,Omegak]=eig(Cov_Yk);
    Yk_Train=Y_set_Tr((k-1)*tau_tol+1:k*tau_tol,:);
    parfor n=1:tau_tol
        Yn_Train=Yk_Train(n,:);
        if bit_allo(n,k)~=0
            [codebook_Y{k,n,:},~]=VQ_Lloyd(Yn_Train,bit_allo(n,k));    
        end
    end
    % N_com=4; % num of time samples sharing a common codebook
    % n_cb=floor(tau_tol/N_com);
    % res_t=tau_tol-n_cb*N_com;
    % for t=1:n_cb
    %     Yn_Train=reshape(Yk_Train((t-1)*N_com+1:t*N_com,:),1,[]);
    %     if bit_allo((t-1)*N_com+1,k)~=0
    %         [codebook_Y{k,(t-1)*N_com+1,:},~]=VQ_Lloyd(Yn_Train,bit_allo((t-1)*N_com+1,k));    
    %         for i=2:N_com
    %             codebook_Y{k,(t-1)*N_com+i,:}=codebook_Y{k,(t-1)*N_com+1,:};
    %         end
    %     end
    % end
    % Yn_Train=reshape(Yk_Train(n_cb*N_com+1:end,:),1,[]);
    % if res_t~=0
    %     if bit_allo(n_cb*N_com+1,k)~=0
    %         [codebook_Y{k,n_cb*N_com+1,:},~]=VQ_Lloyd(Yn_Train,bit_allo(n_cb*N_com+1,k));    
    %         for i=2:res_t
    %             codebook_Y{k,n_cb*N_com+i,:}=codebook_Y{k,n_cb*N_com+1,:};
    %         end
    %     end
    % end
end

%Quantization
Alpha_quan_set=zeros(N*K,num_test);
Y2_quan_set=zeros(tau2*K,num_test);
Y2_set=zeros(tau2*K,num_test);
MSE_Y_quan=zeros(K,1);
pow_Y=zeros(K,1);
NMSE_Y_quan=zeros(K,1);
for k=1:K
    Yk_Test=Y_set_Test((k-1)*tau_tol+1:k*tau_tol,:);
    Yk_quan_w=zeros(tau_tol,num_test);
    for n=1:tau_tol
        % if n<=N
        %     Yn_seg=Yk_Test(n,:)./sqrt(scale_fac);
        %     if ~isempty(codebook_Y{k,n,:})
        %         Yn_quan=Quantizer(Yn_seg,codebook_Y{k,n,:});
        %         Yk_quan_w(n,:)=Yn_quan.*sqrt(scale_fac);
        %     end
        % else
            Yn_seg=Yk_Test(n,:);
            if ~isempty(codebook_Y{k,n,:})
                Yn_quan=Quantizer(Yn_seg,codebook_Y{k,n,:});
                Yk_quan_w(n,:)=Yn_quan;
            end
        % end

    end
    Yk_quan=Yk_quan_w;
    Alpha_quan_set((k-1)*N+1:k*N,:)=Yk_quan(1:N,:);
    Y2_quan_set(k:K:K*(tau2-1)+k,:)=Yk_quan(N+1:end,:);
    Y2_set(k:K:K*(tau2-1)+k,:)=Yk_Test(N+1:end,:);
    MSE_Y_quan(k,:)=mean(sum(abs(Yk_quan-Yk_Test).^2,1));
    pow_Y(k,:)=mean(sum(abs(Yk_Test).^2,1));  
    NMSE_Y_quan(k,1)=MSE_Y_quan(k,:)./pow_Y(k,:);    
end
NMSE_Y=mean(NMSE_Y_quan);
NMSE_all(1)=NMSE_Y;
MSE_Alpha=sum(abs(Alpha_quan_set-Alpha_set).^2);
pow_Alpha=sum(abs(Alpha_set).^2);
NMSE_Alpha=mean(MSE_Alpha)./mean(pow_Alpha);
NMSE_all(2)=NMSE_Alpha;
MSE_Y2=sum(abs(Y2_quan_set-Y2_set).^2,1);
pow_Y2=sum(abs(Y2_set).^2,1);
NMSE_Y2=mean(MSE_Y2)./mean(pow_Y2);
NMSE_all(3)=NMSE_Y2;   

end

%select new reference channel
% G_ref_set=zeros(M*N,num_test);
% Cov_G_ref=zeros(M*N);
% Gk_set=zeros(M*N*K,num_test);
% coeff_newq_set=zeros(N*K,num_test);
% coeff_new_set=zeros(N*K,num_test);
% Alpha_ref_set=zeros(N,num_test);
% Alpha_refq_set=zeros(N,num_test);
% coeff_es_set=zeros(N*(K-1),num_test);
% coeff_true_set=zeros(N*(K-1),num_test);
% for ii=1:num_test
%     coeff_ii=reshape(coeff_out(:,ii),[N,K-1]);
%     LAMDA=[ones(1,N);(coeff_ii).']; % reflecting channel ratio
%     G1=G1_all((ii-1)*M+1:ii*M,:); %M*N     
%     Gk=repmat(G1,[K,1]).*repelem(LAMDA,M,1); % all users' reflecting channel (KM*N)
%     Gk_set(:,ii)=repmat(G1(:),[K,1]).*[ones(M*N,1);repelem(coeff_out(:,ii),M,1)];
% 
%     Alpha_quan=reshape(Alpha_quan_set(:,ii),[N,K]);
%     Alpha=reshape(Alpha_set(:,ii),[N,K]);
%     ref_idx=zeros(N,1);    
%     G_ref=zeros(M*N,1);
%     coeff_new=zeros(N,K);
%     coeff_new_quan=zeros(N,K);
%     Alpha_ref_quan=zeros(N,1);
%     Alpha_ref=zeros(N,1);
% 
%     for n=1:N
%         [~,ref_idx(n)]=max(abs(Alpha_quan(n,:)));
%         % ref_idx(n)=1;
%         G_ref((n-1)*M+1:n*M,:)=Gk((ref_idx(n)-1)*M+1:ref_idx(n)*M,n);
%         coeff_new(n,:)=Alpha(n,:)./Alpha(n,ref_idx(n));
%         coeff_new_quan(n,:)=Alpha_quan(n,:)./Alpha_quan(n,ref_idx(n));
%         Alpha_ref_quan(n,:)=Alpha_quan(n,ref_idx(n));
%         Alpha_ref(n,:)=Alpha(n,ref_idx(n));
%     end
%     coeff_newq_set(:,ii)=coeff_new_quan(:);
%     coeff_new_set(:,ii)=coeff_new(:);
%     G_ref_set(:,ii)=G_ref;
%     Alpha_ref_set(:,ii)=Alpha_ref;
%     Alpha_refq_set(:,ii)=Alpha_ref_quan;
% 
%     Cov_G_ref=Cov_G_ref+(G_ref_set(:,ii))*(G_ref_set(:,ii))';
% 
%     coeff_es=coeff_newq_set(:,ii);
%     coeff_true=coeff_new_set(:,ii);
%     coeff_es_set(:,ii)=coeff_es(abs(coeff_es)~=1);
%     coeff_true_set(:,ii)=coeff_true(abs(coeff_true)~=1);        
% end
% Cov_G_ref=Cov_G_ref./num_test;
% 
% MSE_coeff=sum(abs(coeff_es_set-coeff_true_set).^2,1);
% pow_coeff=sum(abs(coeff_true_set).^2,1);
% NMSE_coeff=mean(MSE_coeff)./mean(pow_coeff);
% NMSE_all(4)=NMSE_coeff;
% 
% %Estimation    
% Delta_quan_set=zeros(tau2*K+N,num_test);
% Delta_set=zeros(tau2*K+N,num_test);
% Cov_err_Delta=zeros(tau2*K+N,tau2*K+N);
% 
% G_ref_quan=zeros(M*N,num_test);
% Gk_quan_set=zeros(M*N*K,num_test);
% Gk_set_new=zeros(M*N*K,num_test);
% OMEGAq_set=zeros(num_test*(tau2*K+N),M*N);
% for ii=1:num_test
%     LAMDA_quan=(reshape(coeff_newq_set(:,ii),[N,K])).';
% 
%     THETAq=Gen_THETA(K,M,N,tau2,P2,LAMDA_quan,pilot_irs_mode,x2);
%     OMEGAq=[THETAq;kron(eye(N),sqrt(P1)*x1.')]; 
% 
%     Delta_quan_set(:,ii)=[Y2_quan_set(:,ii);Alpha_refq_set(:,ii)];
%     Delta_set(:,ii)=[Y2_pure_set(:,ii);Alpha_ref_set(:,ii)];
%     err_Delta=Delta_quan_set(:,ii)-Delta_set(:,ii);
%     Cov_err_Delta=Cov_err_Delta+err_Delta*err_Delta';
%     OMEGAq_set((ii-1)*(tau2*K+N)+1:ii*(tau2*K+N),:)=OMEGAq;
% end
% Cov_err_Delta=Cov_err_Delta./num_test;
% 
% for ii=1:num_test        
%     OMEGAq=OMEGAq_set((ii-1)*(tau2*K+N)+1:ii*(tau2*K+N),:);
%     G_ref_quan(:,ii)=Cov_G_ref*OMEGAq'/(OMEGAq*Cov_G_ref*OMEGAq'+Cov_err_Delta)*Delta_quan_set(:,ii);
%     % G_ref_quan(:,ii)=OMEGAq\Delta_quan_set(:,ii);
%     Gk_quan_set(:,ii)=repmat(G_ref_quan(:,ii),[K,1]).*repelem(coeff_newq_set(:,ii),M,1);
%     Gk_set_new(:,ii)=repmat(G_ref_set(:,ii),[K,1]).*repelem(coeff_new_set(:,ii),M,1);
% end    
% MSE_G_ref=sum(abs(G_ref_quan-G_ref_set).^2,1);
% pow_G_ref=sum(abs(G_ref_set).^2,1);
% MSE_Gk=sum(abs(Gk_quan_set-Gk_set_new).^2,1);
% pow_Gk=sum(abs(Gk_set_new).^2,1);
% NMSE_G_ref=mean(MSE_G_ref)./mean(pow_G_ref);
% NMSE_all(5)=NMSE_G_ref;
% NMSE_Gk=mean(MSE_Gk)./mean(pow_Gk);
% NMSE_all(6)=NMSE_Gk;

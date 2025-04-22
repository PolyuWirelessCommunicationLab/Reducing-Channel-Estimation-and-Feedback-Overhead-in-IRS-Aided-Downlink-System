function THETA=Gen_THETA(K,M,N,tau2,P2,LAMDA,mode,x2)

F=dftmtx(max(N,tau2));
switch mode
    case 1 %cannot work
        [Phi2,A2]=pdc_mtx2(K,M,N,tau2);
        PHI=zeros(tau2,N);
        X=zeros(tau2,M-1);
        for t=1:tau2
            phi2=Phi2(t,:);
            PHI(t,phi2(1):phi2(end))=1;
        end
        X(1:A2(1),1)=1;
        for m=2:M-1
            X(1+sum(A2(1:m-1)):sum(A2(1:m)),m)=1;
        end
        % X=[X,zeros(tau2,1)];
        THETA_bar_s=zeros(tau2*K,(M-1)*N);
        Omegat=zeros(tau2*K,N);
        for t=1:tau2
            IRS_ON=repmat(PHI(t,:),K,1);   
            Omegat((t-1)*K+1:t*K,:)=LAMDA.*IRS_ON;
            THETA_bar_s((t-1)*K+1:t*K,:)=kron(Omegat((t-1)*K+1:t*K,:),sqrt(P2)*X(t,:));
        end
        THETA_bar=[THETA_bar_s,zeros(tau2*K,N)];
        THETA=zeros(tau2*K,M*N);
        for m=1:M
           for n=1:N
               THETA(:,(n-1)*M+m)=THETA_bar(:,(m-1)*N+n);
           end
        end
    case 2
        % THETA_min=prodmtx(N,M,K,LAMDA,P2);
        THETA_min=[];
        tau2_min=0;%ceil(N*(M-1)/K);
        tau2_res=tau2-tau2_min;
        THETA_res=[];
        if tau2_res>0                
            Phi2=F(tau2_min+1:end,1:N);
            x2_res=x2(tau2_min+1:end,:);
            % Phi2=randi([0 1],[tau2_res,N]);
            % x2_res=randi([0 1],[tau2_res,M]);
            for t=1:tau2_res
                THETA_t=zeros(K,M*N);
                for k=1:K
                    THETA_t(k,:)=kron(LAMDA(k,:).*Phi2(t,:),sqrt(P2)*x2_res(t,:));
                end
                THETA_res=[THETA_res;THETA_t];
            end
        end
        THETA=[THETA_min;THETA_res];
    case 3        
        Phi2=F(1:tau2,1:N);
        THETA=zeros(tau2*K,M*N);
        for t=1:tau2
            THETA_t=zeros(K,M*N);
            for k=1:K
                THETA_t(k,:)=kron(LAMDA(k,:).*Phi2(t,:),sqrt(P2)*x2(t,:));
            end
            THETA((t-1)*K+1:t*K,:)=THETA_t;
        end 
end


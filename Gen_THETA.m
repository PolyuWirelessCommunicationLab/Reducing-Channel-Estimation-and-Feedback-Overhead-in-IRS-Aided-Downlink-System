function THETA=Gen_THETA(K,M,N,tau2,P2,LAMDA,x2)

F=dftmtx(max(N,tau2))/sqrt(max(N,tau2));
THETA_min=[];
tau2_min=0;%ceil(N*(M-1)/K);
tau2_res=tau2-tau2_min;
THETA_res=[];
if tau2_res>0                
    Phi2=F(tau2_min+1:end,1:N);
    x2_res=x2(tau2_min+1:end,:);
    for t=1:tau2_res
        THETA_t=zeros(K,M*N);
        for k=1:K
            THETA_t(k,:)=kron(LAMDA(k,:).*Phi2(t,:),sqrt(P2)*x2_res(t,:));
        end
        THETA_res=[THETA_res;THETA_t];
    end
end
THETA=[THETA_min;THETA_res];



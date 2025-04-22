function bit_allo_op=bit_allo_allrev(tau_tol, tau1, K, Cov_Y, Bits_tol)
% bit allication optimization
quan_noise_pow=zeros(tau_tol,K);
bit_allo_op=zeros(tau_tol,K);
max_iter=100;
for k=1:K
    Cov_Yk=Cov_Y((k-1)*tau_tol+1:k*tau_tol,:);
    pow_Yk=diag(Cov_Yk);
    sigma_q=1e-2*rand(tau_tol,1);
    %SCA optimize
    distortion_r=zeros(1,max_iter);
    for iter=1:max_iter
        rate_diff=-1e2;
        %bisection
        lambda=[1e-20 1e-6];
        t=1;
        while abs(rate_diff)>1e-9 
            lambda_mid=mean(lambda);
            delta=pow_Yk.^2+4.*lambda_mid.*(sigma_q+pow_Yk).^2./(log(2)*pow_Yk);
            sigma_q_new=real(1/2*(-pow_Yk+sqrt(delta)));
            rate_diff=sum(log2(1+pow_Yk./sigma_q_new))-Bits_tol;
            if rate_diff>0
                lambda(1)=lambda_mid;
            else
                lambda(2)=lambda_mid;
            end 
            t=t+1;
            if t>1000
                break;
            end
        end     
        sigma_q=sigma_q_new;
        if iter>1&&distortion_r(iter-1)<=sum(pow_Yk)-sum(pow_Yk.^2./(pow_Yk+sigma_q))
            break;
        else
            distortion_r(iter)=sum(pow_Yk)-sum(pow_Yk.^2./(pow_Yk+sigma_q));       
        end
    end
    bit_allo=round(log2(1+pow_Yk./sigma_q));
    bit_allo_op(:,k)=bit_allo;
    quan_noise_pow(:,k)=sigma_q;
end
end

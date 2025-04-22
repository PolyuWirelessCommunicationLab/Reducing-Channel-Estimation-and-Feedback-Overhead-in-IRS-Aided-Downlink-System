function bit_allo_op=bit_allo_twosteps(tau_tol, tau1, K, Cov_Y, Bits_tol)
quan_noise_pow=zeros(tau_tol,K);
bit_allo_op=zeros(tau_tol,K);
max_iter=100;
for k=1:K
    Cov_Yk=Cov_Y((k-1)*tau_tol+1:k*tau_tol,:);
    [Pk,~]=eig(Cov_Yk);
    pow_Yk_all=diag(Cov_Yk);
    Bits(1)=round(0.5*Bits_tol);
    Bits(2)=Bits_tol-Bits(1);
    for step=1:2
        if step==1
            tau=tau1;
            pow_Yk=pow_Yk_all(1:tau1);
            endid=tau1;
        else
            tau=tau_tol-tau1;
            pow_Yk=pow_Yk_all(tau1+1:end);
            endid=tau_tol;
        end
        sigma_q=1e-2*rand(tau,1);
        %SCA optimize
        distortion_r=zeros(1,max_iter);
        for iter=1:max_iter
            % cvx_begin
            %     variable sigma_q_new(tau)
            %     minimize (sum(pow_Yk)-sum(pow_Yk.^2./(pow_Yk+sigma_q)) + (pow_Yk.^2./(sigma_q+pow_Yk)).'*(sigma_q_new-sigma_q) )
            %     subject to
            %         sum((rel_entr(1+sigma_q_new./pow_Yk,sigma_q_new./pow_Yk)+rel_entr(sigma_q_new./pow_Yk,1+sigma_q_new./pow_Yk))./log(2)) <= Bits(step);
            % cvx_end
    
            rate_diff=-1e2;
            %bisection
            lambda=[1e-20 1e-6];
            t=1;
            while abs(rate_diff)>1e-9 
                lambda_mid=mean(lambda);
                delta=pow_Yk.^2+4.*lambda_mid.*(sigma_q+pow_Yk).^2./(log(2)*pow_Yk);
                sigma_q_new=real(1/2*(-pow_Yk+sqrt(delta)));
                rate_diff=sum(log2(1+pow_Yk./sigma_q_new))-Bits(step);
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
        bit_allo=floor(log2(1+pow_Yk./sigma_q));
        ii=1;
        [~,sortid]=sort(pow_Yk,'descend');
        while sum(bit_allo)<Bits(step)       
            bit_allo(sortid(ii))=bit_allo(sortid(ii))+1;
            ii=ii+1;
        end
        % bit_allo_op(:,k)=bit_allo;
        bit_allo_op((step-1)*tau1+1:endid,k)=bit_allo;
        quan_noise_pow((step-1)*tau1+1:endid,k)=sigma_q;
    end    
end
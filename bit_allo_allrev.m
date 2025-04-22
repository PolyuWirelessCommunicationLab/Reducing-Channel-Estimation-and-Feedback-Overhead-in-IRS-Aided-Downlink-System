function bit_allo_op=bit_allo_allrev(tau_tol, tau1, K, Cov_Y, Bits_tol)
quan_noise_pow=zeros(tau_tol,K);
bit_allo_op=zeros(tau_tol,K);
max_iter=100;
for k=1:K
    Cov_Yk=Cov_Y((k-1)*tau_tol+1:k*tau_tol,:);
    % [Pk,~]=eig(Cov_Yk);
    pow_Yk=diag(Cov_Yk);
    scale_fac=mean(pow_Yk(1:tau1))/mean(pow_Yk(tau1+1:tau_tol));
    % pow_Yk(1:tau1)=pow_Yk(1:tau1).*(scale_fac./(P1./P2));
    sigma_q=1e-2*rand(tau_tol,1);
    %SCA optimize
    distortion_r=zeros(1,max_iter);
    for iter=1:max_iter
        % cvx_begin
        %     variable sigma_q_new(tau_tol)
        %     minimize (sum(pow_Yk)-sum(pow_Yk.^2./(pow_Yk+sigma_q)) + (pow_Yk.^2./(sigma_q+pow_Yk)).'*(sigma_q_new-sigma_q) )
        %     subject to
        %         sum((rel_entr(1+sigma_q_new./pow_Yk,sigma_q_new./pow_Yk)+rel_entr(sigma_q_new./pow_Yk,1+sigma_q_new./pow_Yk))./log(2)) <= Bits_tol;
        %         for ii=1:tau_tol
        %             sigma_q_new(ii)<=pow_Yk(ii)./(power(2,3)-1);
        %         end
        % cvx_end

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

        %with minimum bits constraint
        % min_Bits=2;
        % unknowns0=1e-8*rand(2*tau_tol+1,1);
        % fun=@(x)paramfun(x,tau_tol,pow_Yk,sigma_q,Bits_tol,min_Bits);
        % options = optimoptions(@fsolve,'Display','off','Algorithm','trust-region-dogleg');
        % unknowns=fsolve(fun,unknowns0,options);
        % sigma_q_new=unknowns(1:tau_tol);

        sigma_q=sigma_q_new;
        if iter>1&&distortion_r(iter-1)<=sum(pow_Yk)-sum(pow_Yk.^2./(pow_Yk+sigma_q))
            break;
        else
            distortion_r(iter)=sum(pow_Yk)-sum(pow_Yk.^2./(pow_Yk+sigma_q));       
        end
    end
    bit_allo=round(log2(1+pow_Yk./sigma_q));
    % ii=1;
    % [~,sortid]=sort(pow_Yk,'descend');
    % while sum(bit_allo)<Bits_tol       
    %     bit_allo(sortid(ii))=bit_allo(sortid(ii))+1;
    %     ii=ii+1;
    % end
    bit_allo_op(:,k)=bit_allo;
    quan_noise_pow(:,k)=sigma_q;
end
end

function Opt = paramfun(x,tau_tol,pow_Yk,sigma_q,Bits_tol,min_Bits)

    Opt=zeros(2*tau_tol+1,1);
    for i=1:tau_tol
        Opt(i)=(pow_Yk(i).^2./((pow_Yk(i)+sigma_q(i)).^2)+x(tau_tol+1+i)).*x(i).^2+(pow_Yk(i).^3./((pow_Yk(i)+sigma_q(i)).^2)+pow_Yk(i)*x(tau_tol+1+i)).*x(i)-x(tau_tol+1)*pow_Yk(i)./log(2);
    end
    for i=1:tau_tol
        Opt(tau_tol+i)=(x(i)-pow_Yk(i)./(2^min_Bits-1)).*x(tau_tol+1+i);
    end
    Opt(2*tau_tol+1)=sum(log2(1+pow_Yk./x(1:tau_tol)))-Bits_tol;
end
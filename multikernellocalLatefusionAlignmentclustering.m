function [Hstar,WP,gamma,obj] = multikernellocalLatefusionAlignmentclustering(localH,k,lambda,Y,avgH)

num = size(localH, 1); %the number of samples
numker = size(localH, 3); %m represents the number of kernels
maxIter = 100; %the number of iterations
%construct r_p,wp
gamma = ones(numker,1)/(numker);
HP = localH;
WP = zeros(k,k,numker);
opt.disp = 0;
for p=1:numker % m - kernels
    WP(:,:,p) = eye(k);
end

RpHpwp = zeros(num,k); % k - clusters, N - samples


flag = 1;
iter = 0;
while flag
    for p=1:numker
    % Rp2Hpwp = Rp2Hpwp + gamma(p)^2*Hp(:,:,p)*wp(:,:,p);
        RpHpwp = RpHpwp + gamma(p)*(HP(:,:,p)*WP(:,:,p));
    end
    iter = iter +1;
    %the first step-- optimize H_star, given Rp,wp
    %equivalent to maximize(sum(rp2*tr(H_star'*Hp*wp)))
    %SVD of Hp'*H_star
    UU = RpHpwp + lambda * avgH;
    [Uh,Sh,Vh] = svd(UU,'econ');
    Hstar = Uh*Vh';
%     Hstar = zeros(num,k);
%     for ipk =1:num
%         [val, indx] = max(UU(ipk,:));
%         Hstar(ipk,indx) = 1;
%                                                                                                                                                                                                                                                                                                                                                                                                                                                                    
%     end
    %the second step-- optimizing wp, given Rp,H_star equivalent to maximize(sum(2*tr(wp'*Hp'*H_star)))
    %SVD of Hp'*H_star, we will get one wp for each for_loop as follows
    for p=1:numker
        if gamma(p)>1e-4
            TP = gamma(p)*HP(:,:,p)'*(Hstar);
%             WP0 = zeros(k);
%             for ipk =1:k
%                 [val, indx] = max(TP(:,ipk));
%                 e00 = zeros(k,1);
%                 e00(indx) =1;
%                 WP0(:,ipk) = e00;
%             end
%             WP(:,:,p) = WP0;
            [Up,Sp,Vp] = svd(TP,'econ');
            WP(:,:,p) = Up*Vp';
        end  
    end  
    
    
    coef = zeros(1,numker);
    
    for p=1:numker
        coef(1,p) = trace((Hstar)'*HP(:,:,p)* WP(:,:,p)); 
    end
    
    gamma = coef/norm(coef,2);
%%    gamma = ones(numker,1)/(numker);
    
    RpHpwpnew = zeros(num,k);
    for p=1:numker
        RpHpwpnew = RpHpwpnew + gamma(p)*(HP(:,:,p)*WP(:,:,p));
    end
    

    obj(iter) = trace(Hstar'*RpHpwpnew + lambda * Hstar'*avgH);
    
%     res9 = myNMIACC(Hstar,Y,k);
%     acc(iter) = res9(1);
%     nmi(iter)= res9(2);
%     pur(iter) = res9(3);
    
    if (iter>2) && (abs((obj(iter-1)-obj(iter))/(obj(iter-1)))<1e-6 || iter>maxIter)
%      if iter==maxIter
       flag =0;
    end
end
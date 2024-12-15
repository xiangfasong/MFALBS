function [idx,obj] = MFALBS(X,W,G,B,E,H,F,param)
%% Problem
% min_W,S,B,E=||W'X-BE'||_sigma+lambda*||W||_21+alhpa*Tr(W'XLX'W)+beta*||S-G||^F_2
% s.t. B'B=I, E'E = I, E>=0
alpha = param.alpha;
lambda = param.lambda;
beta = param.beta;
sigma = param.sigma;
Y1 = zeros (size(W,2),size(X,2));%m*n
Y2 = zeros (size(X,2),size(B,2));%n*c
rho = 1e2;                              
NITER = 50;
S = (G+G')/2;
L_D = diag(sum(S));
L = L_D-S;
[dimX,nX] = size(X);
Nsap = nX;
IIT = ones(Nsap,Nsap);
ac1 = ones(Nsap,1);
P = ones(dimX,1);
Q = ones(size(F,1),1);
P = spdiags(P,0,dimX,dimX);
Q = spdiags(Q,0,size(F,1),size(F,1));
if (dimX > nX) 
    for iter=1:NITER
        %Update B
        Temp1 = E'*(X'*W-F'-Y1'./rho);
        [a,~,c] = svd(Temp1','econ');
        B = a * c';
        %Update W
        PX = P*X;
        W1 = (2*alpha*L*X'+rho*X')*PX+2*lambda*eye(nX);
        W2 = rho*F+rho*B*E'+Y1;
        W = PX*pinv(W1)*W2';
        % update P
        wc = sqrt(sum(W.*W,2)+eps);
        P = 2.*wc;
        P = spdiags(P,0,dimX,dimX);
        %  Update S 
        Y = X'*W;
        ab=0.5*alpha/beta;
        U = G+ab*(Y*Y');    
        U = (U+U')*0.5;
        UIIT = U*IIT;
        S = U+1/Nsap^2*((Nsap+ac1'*U*ac1)*IIT)-1/Nsap*(UIIT+UIIT');
        S = S-diag(diag(S));
        for in_iter = 1:10
             if min(S(:)) < 0
                 S = max(S,0);           
                 U = S;
                 U = (U+U')*0.5;
                 UIIT = U*IIT;
                 S = U+1/Nsap^2*((Nsap+ac1'*U*ac1)*IIT)-1/Nsap*(UIIT+UIIT');
                 S = S-diag(diag(S));
             else
                 break;
             end
         end
         S = max(S,0);
        S = (S+S')/2;
        L_D = diag(sum(S));
        L = L_D-S;
        % Update E
        Temp2 = (W'*X-F-Y1./rho)'*B+H+Y2./rho;
        [a,~,c] = svd(Temp2,'econ');
        E = a * c';
        % Update H
        H = max(0,E-Y2./rho);
        % Update F
        Temp3 = W'*X-B*E'-Y1./rho;
        F = pinv(2*Q + rho * eye(size(F,1)))*(rho*Temp3);
        %Update Q
        f = sqrt(sum(F.*F,2));
        Q = 0.5*(1+sigma)*(f+2*sigma)./((f+sigma).^2);
        Q = spdiags(Q,0,size(F,1),size(F,1));
        %Update Y1,Y2,rho
        Y1 = Y1 + rho * (F - W'*X + B * E');
        Y2 = Y2 + rho * (H - E);
%         rho = min(rho * pho_rho, max_rho);
        ec= W'*X-B*E';
        sigmaec=sum(((1+sigma).*sum(ec.*ec,2))./(sqrt(sum(ec.*ec,2))+sigma));
        obj(iter) = sigmaec + lambda * norm_21(W) + alpha * trace(W'*X*L*X'*W) + beta * norm(S-G,'fro')^2;
    end
else
    for iter=1:NITER
        %Update B
        Temp1 = E'*(X'*W-F'-Y1'./rho);
        [a,~,c] = svd(Temp1','econ');
        B = a * c';
        %Update W
        W1 = 2*alpha*X*L*X'+rho*(X*X')+2*lambda*P;
        W2 = rho*F+rho*B*E'+Y1;
        W = pinv(W1)*X*W2';
        % update P
        wc = sqrt(sum(W.*W,2)+eps);
        P = 0.5./wc;
        P = spdiags(P,0,dimX,dimX);
        % Update S 
        Y = X'*W;
        ab=0.5*alpha/beta;
        U = G+ab*(Y*Y');    
        U = (U+U')*0.5;
        UIIT = U*IIT;
        S = U+1/Nsap^2*((Nsap+ac1'*U*ac1)*IIT)-1/Nsap*(UIIT+UIIT');
        S = S-diag(diag(S));
        for in_iter = 1:10
             if min(S(:)) < 0
                 S = max(S,0);           
                 U = S;
                 U = (U+U')*0.5;
                 UIIT = U*IIT;
                 S = U+1/Nsap^2*((Nsap+ac1'*U*ac1)*IIT)-1/Nsap*(UIIT+UIIT');
                 S = S-diag(diag(S));
             else
                 break;
             end
         end
        S = max(S,0);
        S = (S+S')/2;
        L_D = diag(sum(S));
        L = L_D-S;
        % Update E
        Temp2 = (W'*X-F-Y1./rho)'*B+H+Y2./rho;
        [a,~,c] = svd(Temp2,'econ');%Temp2的奇异的左特征向量和右特征向量就是Temp2'的奇异的右特征向量和左特征向量
        E = a * c';
        % 更新H
        H = max(0,E-Y2./rho);
        % Update F
        Temp3 = W'*X-B*E'-Y1./rho;
        F = pinv(2*Q + rho * eye(size(F,1)))*(rho*Temp3);
        %Update Q
        f = sqrt(sum(F.*F,2));
        Q = 0.5*(1+sigma)*(f+2*sigma)./((f+sigma).^2);
        Q = spdiags(Q,0,size(F,1),size(F,1));
        %Update Y1,Y2,rho
        Y1 = Y1 + rho * (F - W'*X + B * E');
        Y2 = Y2 + rho * (H - E);
%         rho = min(rho * pho_rho, max_rho);
        ec= W'*X-B*E';
        sigmaec=sum(((1+sigma).*sum(ec.*ec,2))./(sqrt(sum(ec.*ec,2))+sigma));
        obj(iter) = sigmaec + lambda * norm_21(W) + alpha * trace(W'*X*L*X'*W) + beta * norm(S-G,'fro')^2;
    end
end
[~,idx] = sort(sum(W.*W,2),'descend');
end


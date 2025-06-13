%% Get prior estimates A0 and B0
%D0 uncertainty quantification
% mul=1;
% cbar=(sigma_w^2)*chi2inv(delta,(nx*(nx+nu)));
% D0=(mul*1e5)*eye(nx+nu);
% d0=(cbar/sigma_w^2)*D0;
% kprior=kron(d0,eye(nx));
% covprior=inv(kprior);

count=1;
while count>0
    theta0=(mvnrnd(thetatr,kprior))';
    tr=(thetatr-theta0)'*D0*(thetatr-theta0);
    if tr<1
%     if abs(tr-cbar)<1e-2
        count=0;
    end
end

p0=reshape(theta0,[nx,nphi]);
A0=p0(:,1:nx);
B0=p0(:,nx+1);

%% Compute \hat{V}, \hat{Y}, ||\hat{V}-V_tr||, ||Y_tr||

Vh=[]; Vtr=[]; Ytr=[]; Yh=[];

for i=1:nomega
    vt=inv(exp(1i*2*pi*freqs(i))*eye(nx)-A0)*B0;
    vi=[vt;eye(nu)];
    Vh=[Vh,vi];
end
nVh=norm(Vh); % || V_hat ||

for i=1:nomega
    vt=inv(exp(1i*2*pi*freqs(i))*eye(nx)-A)*B;
    vi=[vt;eye(nu)];
    Vtr=[Vtr,vi];
end
nVtr=norm(Vtr); % || V_tr ||

for i=1:nomega
   yt=inv(exp(1i*2*pi*freqs(i))*eye(nx)-A);
   yi=[yt;zeros(nu,nx)];
   Ytr=[Ytr,yi];
end
nYtr=norm(Ytr); % || Y_tr ||

for i=1:nomega
   yt=inv(exp(1i*2*pi*freqs(i))*eye(nx)-A0);
   yi=[yt;zeros(nu,nx)];
   Yh=[Yh,yi]; 
end
nYh=norm(Yh); % || Y_h ||

Vtilde_tr=Vh-Vtr;
nVtilde_tr=norm(Vtilde_tr); % || V_hat - V_tr ||


%% true H2 norm - closed loop system

% cvx_begin SDP quiet
% variable N(nx,nx) symmetric
% variable M(nu,nx)
% variable Z(nx,nx) symmetric
% % variable Z(1)
% variable g(1)
% S2=[A*N+B*M, eye(nx);zeros(nx), zeros(nx)];
% S3=blkdiag(-N,-eye(nx));
% S1=blkdiag(-N,-(g)*eye(nx));
% S=[S1,S2';S2,S3];
% minimize g
% subject to
% S<=0;
% [N N'*C';C*N Z]>=0;
% trace(Z)<=g;
% % Z<=g;
% g>=0;
% cvx_end


%% 
% C=eye(nx);
% cvx_begin SDP quiet
% variable N(nx,nx) symmetric
% variable M(nu,nx)
% variable Z(nx,nx) symmetric
% % variable Z(1,1) symmetric
% variable g0(1)
% lambda=20;
% S2=[A0*N+B0*M, eye(nx), eye(nx);[N;M], zeros(nphi,nx), zeros(nphi,nx);zeros(nx), zeros(nx), zeros(nx)];
% S3=blkdiag(-N,-(1/lambda)*(D0),-eye(nx));
% S1=blkdiag(-N,-lambda*eye(nx),-(g0)*eye(nx));
% S=[S1,S2';S2,S3];
% 
% minimize g0
% subject to
% S<=0;
% % trace(N)<=g0;
% [N N'*C';C*N Z]>=0;
% trace(Z)<=g0;
% g0>=0;
% 
% cvx_end


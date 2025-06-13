%% Sample based bound for gamma_v, gamma_y and l_1

theta0=[A0(:);B0];
thetatr=[A(:);B];
thetas=[];
params=[];
n=0;
Ns=500;

while n<Ns
    t1=(mvnrnd(theta0,kprior))';
    t2=(theta0-t1)'*D0*(theta0-t1);
    if t2<=1
        thetas=[thetas,t1];
        paramt=reshape(t1,[nx,nphi]);
        params=[params;paramt];
        n=n+1;
    end
end


%% Determine ||V_tilde||, ||Y_tr|| from samples

Yhs=[];
Vhs=[];
NY=[]; NV=[]; NVtilde=[];
for i=1:nx:nx*n
    pa=params(i:i+nx-1,1:nx);
    pb=params(i:i+nx-1,nx+1:nphi);
    for j=1:nomega
        ytemp=inv(exp(1i*2*pi*freqs(j))*eye(nx)-pa);
        vtemp=inv(exp(1i*2*pi*freqs(j))*eye(nx)-pa)*pb;
        yi=[ytemp;zeros(nu,nx)];
        vi=[vtemp;eye(nu)];
        Yhs=[Yhs,yi];
        Vhs=[Vhs,vi];
    end
    NV=[NV,norm(Vh-Vhs)];
    NVtilde=[NVtilde,norm((Vh-Vhs)*((Vh'*Vh)^(-1/2)))];
    NY=[NY,norm(Yhs)];
    Yhs=[]; Vhs=[];
end
% nYhs=mean(NY); %this value is gamma_y
% nVhs=mean(NV);
% nVtilde=mean(NVtilde); %this value is gamma_v
gammav=max(NV);
gammay=max(NY);
gammav1=max(NVtilde);
% gammav1=gammav*norm(inv(sqrt(Vh*Vh')));

%% Determine l_1=||W||

sigma=sigma_w/(sqrt(T));
l1d=sigma*sqrt(chi2inv(1-delta,nx));
% l12d=l1d*l1d;
l=l1d*gammay;


%% Determine \Gamma_v

Vhs=[];
cvx_begin SDP quiet
% cvx_solver mosek
cvx_precision high

variable Gammav((nx+nu),(nx+nu)) semidefinite

minimize trace(Gammav)
subject to
for i=1:nx:nx*n
    pa=params(i:i+nx-1,1:nx);
    pb=params(i:i+nx-1,nx+1:nphi);
    for j=1:nomega
        vtemp=inv(exp(1i*2*pi*freqs(j))*eye(nx)-pa)*pb;
        vi=[vtemp;eye(nu)];
        Vhs=[Vhs,vi];
    end
    Gammav-((Vh-Vhs)*(Vh-Vhs)')>=0;
    Vhs=[];
end

cvx_end


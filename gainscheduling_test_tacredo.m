%% RUN this after prior_new_tacredo
%%
cvx_begin SDP quiet
variable N(nx,nx) symmetric
variable M(nu,nx)
variable Ks(nu,nx)
variable Dt(nx+nu,nx+nu) semidefinite
variable Z(nx,nx) symmetric
variable gammap(1)
minimize gammap
subject to
lambda=17;
lambdas=100000;
% gammap=3.05;
Sgs_1=blkdiag(-N,-lambda*eye(nx),-lambdas*eye(nx),-(gammap)*eye(nx));
Sgs_3=blkdiag(-N,-(1/lambda)*D0,-(1/lambdas)*(D0+Dt),-1*eye(nx));
Sgs_2=[A0*N+B0*M, eye(nx)+B0*Ks, eye(nx), eye(nx);...
    [N;M], [zeros(nx,nx);Ks], zeros(nphi,nx), zeros(nphi,nx);...
    [N;M], [zeros(nx,nx);Ks], zeros(nphi,nx), zeros(nphi,nx);...
    zeros(nx), zeros(nx), zeros(nx), zeros(nx)];
% last row in Sgs_2 eye(nx)*N = C*N

Sgs=[Sgs_1, Sgs_2'; Sgs_2, Sgs_3];

Sgs<=0;
[N N'*C';C*N Z]>=0;
trace(Z)<=gammap;
% trace(N)<gammap;
cvx_end

% disp(trace(Dt));

%%

cvx_begin SDP quiet
variable N(nx,nx) symmetric
variable M(nu,nx)
variable Ks(nu,nx)
variable Z(nx,nx) symmetric
% variable Dt(nx+nu,nx+nu) semidefinite
variable g00(1)
minimize g00
subject to
% gammap=50;
lambda=17;
Sgs_1=blkdiag(-N,-lambda*eye(nx),-(g00)*eye(nx));
Sgs_3=blkdiag(-N,-(1/lambda)*D0,-1*eye(nx));
Sgs_2=[A0*N+B0*M, eye(nx)+B0*Ks, eye(nx);...
    [N;M], [zeros(nx,nx);Ks], zeros(nphi,nx);...
    zeros(nx), zeros(nx), zeros(nx)];
% last row in Sgs_2 eye(nx)*N = C*N

Sgs=[Sgs_1, Sgs_2'; Sgs_2, Sgs_3];

Sgs<=0;
[N N'*C';C*N Z]>=0;
trace(Z)<=g00;
% trace(N)<gammap;
cvx_end
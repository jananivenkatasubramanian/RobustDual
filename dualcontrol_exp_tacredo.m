%% Solve joint problem of exploration and gain-scheduled control

cvx_begin SDP quiet
% cvx_solver mosek
% cvx_precision high

%% VARIABLES
%Dt lower bound
% variable Dt((nx+nu),(nx+nu)) semidefinite
% Dt == semidefinite(nphi);

variable u(nu,nomega)

Ue=[];
for i=1:nomega
    Ue=blkdiag(Ue,u(:,i));
end

variable gammae(1)
variable tau(1)
%---------------------------------------------
% Exploration energy

Se_energy=[gammae, ones(1,nomega)*(Ue'); Ue*ones(nomega,1), gammae*eye(nomega)];

%---------------------------------------------
%Exploration LMI

Pe11=(1-eps)*(Ue*L'+L*Ue'-L*L');
Pe22=-((1-eps)/eps)*(l^2)*eye(nphi) - (nomega/T)*Dt;
Pe=blkdiag(Pe11,Pe22);

% Qe=[-eye(nphi), Vh';Vh, (gammav^2)*eye(nphi) - Vh*Vh'];
Qe=[-eye(nomega), Vh';Vh, Gammav - Vh*Vh'];
% Qe11=eye(nphi);
% Qe22=-(1-2*gammav1)*(Vh*Vh');
% Qe=blkdiag(Qe11,Qe22);

Se=Pe-tau*Qe;
%---------------------------------------------
% Gain Scheduling LMI
% Sgs_1=blkdiag(-N,-lambda*eye(nx),-lambdas*eye(nx),-(gammap)*eye(nx));
% Sgs_3=blkdiag(-N,-(1/lambda)*D0,-(1/lambdas)*(D0+Dt),-1*eye(nx));
% Sgs_2=[A0*N+B0*M, eye(nx)+B0*Ks, eye(nx), eye(nx);...
%     [N;M], [zeros(nx,nx);Ks], zeros(nphi,nx), zeros(nphi,nx);...
%     [N;M], [zeros(nx,nx);Ks], zeros(nphi,nx), zeros(nphi,nx);...
%     zeros(nx), zeros(nx), zeros(nx), zeros(nx)];
% 
% Sgs=[Sgs_1, Sgs_2'; Sgs_2, Sgs_3];

%---------------------------------------------
minimize gammae
subject to
tau>=0;
% Dt(1,1)-1e6>=0;
Se_energy>=0;
Se>=0;
% Sgs<=0;
% [N N'*C';C*N Z]>=0;
% trace(Z)<=gammap;


cvx_end


% Ut=Ue;
% L=V*(1)*Ue;
L=Ue;
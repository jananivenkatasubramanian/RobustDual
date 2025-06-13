%% Initialization - Sequential learning and control
clear all;
clc;
% system order
% nx=3;
nx=4;
nu=1;
% nphi=nx+nu;
nphi=5;
nomega=10;

% probability level
delta=0.01;
eps=0.5;
% sys=drss(nx,nx,nu);
% A=sys.A;
% B=sys.B;
% C=sys.C;
% D=sys.D;
% h=hinfnorm(sys);

% A=[0.49 0.49 0; 0 0.49 0.49; 0 0 0.49];
% B=[1; 0; 0.49];  
A=[0.49 0.49 0 0;0 0.49 0.49 0; 0 0 0.49 0.49;0 0 0 0.49];
B=[0; 0; 0; 0.49];
% C=ones(1,nx);
C=eye(nx);
% C=blkdiag(1,0,0,0);
D=zeros(nx,nu);
% PP=ss(A,B,C,D,1);
thetatr=[A(:);B];
% h=hinfnorm([A,B,C,D]);

% h=h(1);

% u and w parameters to obtain prior
mu_u=0;
sigma_u=sqrt(2000);

mu_w=0;
sigma_w=1;

% number of trials, number of rollouts
N_t=1;
N_r=1000;


% Exploration time
T=100;

% Omega_T set of possible frequencies
Omega_T=[];
for i=1:T
    Omega_T(i)=(i-1)/T;
end

% freqs=[];

% Set of nphi selected frequencies
% sf=randperm(T,nphi);
% for i=1:nphi
%     freqs(i)=Omega_T(sf(i)); %randomly selected frequencies
% end

%uniformly selected frequencies

% freqs=Omega_T(1:T/nphi:T);

% freqs=[0,0.05,0.1,0.2,0.25];

% freqs=[0,0.04,0.08,0.16,0.2];

% freqs=[0.50, 0.53, 0.01, 0, 0.04];

% freqs=[0, 0.1];

% freqs=[0, 0.01, 0.04, 0.51, 0.53];

% freqs=[0,0.05,0.1,0.15,0.2];

freqs=[0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9];

% freqs=[0,0.2,0.4,0.6,0.8];




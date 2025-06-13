%% Random Exploration

X_exp=zeros(nx,T+1);
U_exp=normrnd(0,sigma_u,nu,T);
U_exp=U_exp*norm(Uc)/norm(U_exp);
% U_exp=normrnd(0,3e6,nu,T);
% w_exp=normrnd(mu_w, sigma_w, nx, T);
x0=zeros(nx,1);
X_exp(:,1)=x0;
for m=1:T
    x=A*X_exp(:,m)+B*U_exp(m)+w_exp(:,m);
    X_exp(:,m+1)=x;
end

d=zeros(nphi);
for m=1:T
    d=d+[X_exp(:,m+1);U_exp(:,m)]*[X_exp(:,m+1);U_exp(:,m)]';
end
Dt_rand_exp=d;
% Dt_rand_exp=d/cbar;
% disp("Dt_exp =");
% disp(Dt_exp);
% disp("min eig Dt =");
% disp(min(eig(Dt_exp)));
% 
% disp("max eig Dt =");
% disp(max(eig(Dt_exp)));

ce_re=0;
for m=1:T
    ce_re=ce_re+(U_exp(1,m)^2);
end
% disp("ce_re =");
% disp(ce_re);

p=[mul,min(eig(Dt_rand_exp)),max(eig(Dt_rand_exp)),Dt_rand_exp(1,1),ce_re];



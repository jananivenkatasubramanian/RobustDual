%% Compute targted exploration - control energy
% Dt=Dt100;
% Ue=Ue100_1;
Uc=zeros(nphi,T);
X=zeros(nx,T+1);
tset=1:T;
for i=1:nphi
    cosi=cos(2*pi*freqs(i)*tset);
    Uc(i,:)=Ue(i,i)*cosi;
end
Uc=sum(Uc);

U_ce=zeros(1,T);
for i=1:T
    U_ce(i)=U_ce(i)+(Uc(i)^2);
end
ce_te=sum(U_ce);
% disp("ce_te =");
% disp(sum(U_ce));

% disp("Dt =");
% disp(Dt);

%% verify Dt_targeted_exp


U_te=1*Uc;
x0=zeros(nx,1);
X(:,1)=x0;
for m=1:T
    x=A*X(:,m)+B*U_te(m)+w_exp(:,m);
    X(:,m+1)=x;
end

d2=zeros(nphi);
for m=1:T
    d2=d2+[X(:,m+1);U_te(:,m)]*[X(:,m+1);U_te(:,m)]';
end

% kd1=kron((d2/(sigma_w^2)),eye(nx));
% kd=kron((d0/(sigma_w^2)),eye(nx));
% Dt_texp=d2/cbar;

Dt_texp=d2;

% disp("Dt_targeted_exp =");
% disp(Dt_texp);
% disp("min_eig(Dt_targeted_exp) =");

% disp(min(eig(Dt_texp)));
% disp(max(eig(Dt_texp)));
% disp(Dt_texp>Dt);

pars=[mul,min(eig(Dt_texp)),max(eig(Dt_texp)),Dt_texp(1,1),ce_te];

% pars=[mul,min(eig(Dt_texp(1:4,1:4))),max(eig(Dt_texp(1:4,1:4))),ce_te]

%%


% cost=0;
% cvx_begin quiet
%     variables AT(nx,nx) BT(nx,nu)
%     for k=1:T
%         y=X(:,k+1)-AT*X(:,k)-BT*U_te(:,k);
%         cost=cost+y'*y;
%     end
%     minimize cost
% cvx_end

% -------------------UNCOMMENT FROM HERE-----------------

% cost=0;
% cvx_begin quiet
%     variables thT(nx*nx+nx,1)
%     for k=1:T
%         y1=norm(X(:,k+1)-kron([X(:,k)' U_te(:,k)'],eye(nx))*thT);
%         y2=norm(sqrt(kd)*(thT-theta0));
%         cost=cost+y1+y2;
%     end
%     minimize cost
% cvx_end
% 
% param=reshape(thT,[nx,nx+nu]);
% AT=param(:,1:nx);
% BT=param(:,nx+1);
% 
% thT_0=(thT-theta0)'*(kd1+kd)*(thT-theta0);
% X=[];
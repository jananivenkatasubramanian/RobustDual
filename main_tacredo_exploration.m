% Main - sequential learning and control

run initialize_sequential.m

w_exp=normrnd(mu_w, sigma_w, nx, T);
muldt_targ=[];
muldt_rand=[];
gpge=[];
%%
% for outer=1:5
    % mul=10^(outer-1);
    mul=10^(4);
    for inner=1:1
        cbar=(sigma_w^2)*chi2inv(delta,(nx*nphi));
        cdelta=chi2inv(delta,(nx*nphi));
        D0tilde=(mul*2e2)*eye(nphi);
        D0=kron(D0tilde,eye(nx));
        D0inv=inv(D0);
        kprior=(1/cbar)*D0inv; %simplified since sigma_w=1

        run prior_new_tacredo.m
        run theta_samples_tacredo.m

        % Dt=1e5*eye(nphi);
        Ddes=1e6*eye(nphi);
        Dt=(cbar*Ddes)-D0tilde;
        L=1e4*eye(nomega);
        gtemp=inf;
        for count=1:20
            run dualcontrol_exp_tacredo.m
            if (abs(gtemp-gammae)/gammae)<1e-2
                break;
            end
            % disp(count);
            gtemp=gammae;
        end

        disp('cvx status=');
        disp(cvx_status);
        disp(Ue);
        disp(inner);
        run controlenergy_sequential.m
        muldt_targ=[muldt_targ;pars];
        
        run randexp_sequential.m
        muldt_rand=[muldt_rand;p];
    end
% end

%%
save("slc2_inif.mat","muldt_targ");
save("slc2_inif.mat","muldt_rand","-append");

%% Debug
% Dt=1e5*eye(nphi);
Dt=blkdiag(1e6,0,0,0,0)*eye(nphi);
L=5e4*eye(nphi);
Ue=8e4*eye(nphi);
tau=2e9;
eps=0.5;

%%
L=1e5*eye(nomega);
Ue=1e5*eye(nphi);
tau=1e9;
%%
% Pe11=(1-eps)*(Ue*L'+L*Ue'-L*L');
Pe11=(1-eps)*(Ue*Ue');
Pe22=-((1-eps)/eps)*(l^2)*eye(nphi) - (cbar*nphi/T)*Dt;
Pe=blkdiag(Pe11,Pe22);

% Qe=[-eye(nphi), Vh';Vh, (gammav^2)*eye(nphi) - Vh*Vh'];
Qe=[-eye(nphi), Vh';Vh, Gammav - Vh*Vh'];
% Qe11=eye(nphi);
% Qe22=-(1-2*gammav1)*(Vh*Vh');
% Qe=blkdiag(Qe11,Qe22);

Se=Pe-tau*Qe;

disp(eig(Se));
disp(min(eig(Se)));
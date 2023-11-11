% Code related to Automatica Paper "On the Application of Galerkin Projection based Polynomial Chaos in Linear Systems and Control"
% by LL Evangelisti and H Pfifer
clear all; clc; close all;

%% 
load('projPCECoEx20.mat') % this matrix contains the projected DeltaPi matrix computed by the PolyChaos.jl Julia Toolbox

p = ureal('p',0,'Range',[-1,1]);
A = 0.01*[128*p^2-72*p-32,  295*p^2-199*p+4,   165*p^2-234*p+46; ...
		  -82*p^2-59*p+270, -266*p^2+144*p-73,  -147*p^2-210*p+286; ...
		  70*p^2+296*p-80,  43*p^2+96*p+8,      15*p^2+146*p-251];
%%
maxdeg=16; gvec = zeros(maxdeg+1,1);
for deg = 0:maxdeg
    L = deg + 1;

    gTry = 31; gLow = 0; gUpp = Inf; RelTol = 1e-2; AbsTol = 1e-3;

    while( (gUpp - gLow >= RelTol*gUpp) && (gUpp - gLow >= AbsTol) )
%         Ham_pi = [Api, eye(L*3)/gTry^2; -eye(L*3) -Api'];
        Ham = [A, eye(3)/gTry^2; -eye(3) -A'];
        [M,Delta,BLKSTRUCT] = lftdata(Ham);
        nd = length(Delta);
        M11 = M(1:nd,1:nd);         M12 = M(1:nd,(nd+1):end); 
        M21 = M((nd+1):end,1:nd);   M22 = M((nd+1):end,(nd+1):end);
        DeltaPiRep = kron(eye(BLKSTRUCT(1).Occurrences), DeltaPi(1:L,1:L));
        Ham_pi = lft( DeltaPiRep, kron(M, eye(L)) );
        
        if(min(abs(real(eig(Ham_pi))))<1e-3)
    %     if(sol.x(1) < T)
            gLow = gTry;
            if(isfinite(gUpp))
                gTry = (gLow + gUpp)/2;
            else
                gTry = gTry*2;
            end
        else
            gUpp = gTry;
            gTry = (gLow + gUpp)/2;
        end
    end
    
    gvec(deg+1)=gUpp;
end
%%
figure;
xvec = 1:(maxdeg+1);
plot(xvec,gvec); hold on
plot(1:100,ones(size(1:100)).*46.1216,'r-')
%%
N = 100;
N_mcs = 100;
mwc = 0;
wc_MCSN = zeros(N_mcs,N);
for j = 1:N_mcs
gwc = 0;
fprintf('%j-',j);
if round(j/20)==j/20 fprintf('\n'); end;
for i = 1:N
    delta_s = -1+rand*2;
    As = usubs(A,'p',delta_s); 
    
    gTry = 31; gLow = 0; gUpp = Inf; RelTol = 1e-2; AbsTol = 1e-3;

    while( (gUpp - gLow >= RelTol*gUpp) && (gUpp - gLow >= AbsTol) )
        Ham = [As, eye(3)/gTry^2; -eye(3) -As'];
        
        if(min(abs(real(eig(Ham))))<1e-3)
            gLow = gTry;
            if(isfinite(gUpp))
                gTry = (gLow + gUpp)/2;
            else
                gTry = gTry*2;
            end
        else
            gUpp = gTry;
            gTry = (gLow + gUpp)/2;
        end
    end
    
    gi = gUpp;
    
    if(gi > gwc)
        gwc = gi;
        wc_MCSN(j,i:end) = gwc;
    end
end
end
%%
ylabel('\gamma'); xlabel('Dimension N')
hpwc = plot(1:N,ones(N,1)*max(max(wc_MCSN)),'r--','LineWidth',2);
mean_MCS_wc = sum(wc_MCSN,1)/N_mcs;
hpmean = plot(1:N,mean_MCS_wc,'k-')
std_MCS_wc = sqrt(sum((wc_MCSN-repmat(mean_MCS_wc,N_mcs,1)).^2,1)/(N_mcs-1));
plot(1:N,mean_MCS_wc+std_MCS_wc,'k--');
plot(1:N,mean_MCS_wc-std_MCS_wc,'k--');
%% 
addpath(genpath('Legendre-Gauss-Quadrature-master'))
maxdeg=99; gvec = zeros(maxdeg+1,1);
for deg = 0:maxdeg
    L = deg + 1;
    
    deltavec = legzo(L);
    
    for i=1:length(deltavec)
        
        delta_s = deltavec(i);
        As = usubs(A,'p',delta_s); 

        gTry = 31; gLow = 0; gUpp = Inf; RelTol = 1e-2; AbsTol = 1e-3;

        while( (gUpp - gLow >= RelTol*gUpp) && (gUpp - gLow >= AbsTol) )
            Ham = [As, eye(3)/gTry^2; -eye(3) -As'];

            if(min(abs(real(eig(Ham))))<1e-3)
                gLow = gTry;
                if(isfinite(gUpp))
                    gTry = (gLow + gUpp)/2;
                else
                    gTry = gTry*2;
                end
            else
                gUpp = gTry;
                gTry = (gLow + gUpp)/2;
            end
        end

        gi = gUpp;

        if(gi > gwc)
            gwc = gi;
            wc_MCSN(j,i:end) = gwc;
        end 
    end
    
    gvec(deg+1)=gUpp;
end
%%
plot(1:(maxdeg+1),gvec,'--');
%%
% T = table(gvec,mean_MCS_wc',std_MCS_wc','VariableNames',{'gauss','meanMC','stdMC'});
% writetable(T,'D:\UQ_PCE\CDC23\figures\counterexample_hinf.csv')
%%
n=100; deltavec = linspace(-1,1,n); gvec = zeros(1,n);
for i = 1:n    
    
    delta_s = deltavec(i);
    As = usubs(A,'p',delta_s); 

    gTry = 31; gLow = 0; gUpp = Inf; RelTol = 1e-2; AbsTol = 1e-3;

    while( (gUpp - gLow >= RelTol*gUpp) && (gUpp - gLow >= AbsTol) )
        Ham = [As, eye(3)/gTry^2; -eye(3) -As'];

        if(min(abs(real(eig(Ham))))<1e-3)
            gLow = gTry;
            if(isfinite(gUpp))
                gTry = (gLow + gUpp)/2;
            else
                gTry = gTry*2;
            end
        else
            gUpp = gTry;
            gTry = (gLow + gUpp)/2;
        end
    end

    gi = gUpp;

    if(gi > gwc)
        gwc = gi;
        wc_MCSN(j,i:end) = gwc;
    end 
    
    gvec(i)=gUpp;
end
%%
figure
plot(deltavec,gvec,'-'); hold on
GQP1 = legzo(1);    plot(GQP1,zeros(size(GQP1)),'ksq');     ga1 = interp1(deltavec,gvec,GQP1);      plot(GQP1,ga1,'ksq')
GQP2 = legzo(2);    plot(GQP2,zeros(size(GQP2)),'ko');      ga2 = interp1(deltavec,gvec,GQP2);      plot(GQP2,ga2,'ko')
GQP4 = legzo(4);    plot(GQP4,zeros(size(GQP4)),'k*');      ga4 = interp1(deltavec,gvec,GQP4);      plot(GQP4,ga4,'k*')
GQP8 = legzo(8);    plot(GQP8,zeros(size(GQP8)),'k+');      ga8 = interp1(deltavec,gvec,GQP8);      plot(GQP8,ga8,'k+')
GQP16 = legzo(16);  plot(GQP16,zeros(size(GQP16)),'kx');    ga16 = interp1(deltavec,gvec,GQP16);    plot(GQP16,ga16,'kx')
ylabel('\gamma'); xlabel('\delta')
% GQP = legzo(4); plot(GQP,zeros(size(GQP)),'x')
%%
% T = table(deltavec',gvec', 'VariableNames',{'d','g'});
% writetable(T,'D:\UQ_PCE\Diss\figures\counterexample_hinf_sample.csv')
% 
% T = table(GQP1,ga1, 'VariableNames',{'d1','g1'});
% writetable(T,'D:\UQ_PCE\Diss\figures\counterexample_hinf_sample1.csv')
% 
% T = table(GQP2',ga2', 'VariableNames',{'d2','g2'});
% writetable(T,'D:\UQ_PCE\Diss\figures\counterexample_hinf_sample2.csv')
% 
% T = table(GQP4',ga4', 'VariableNames',{'d4','g4'});
% writetable(T,'D:\UQ_PCE\Diss\figures\counterexample_hinf_sample4.csv')
% 
% T = table(GQP8',ga8', 'VariableNames',{'d8','g8'});
% writetable(T,'D:\UQ_PCE\Diss\figures\counterexample_hinf_sample8.csv')
% 
% T = table(GQP16',ga16', 'VariableNames',{'d16','g16'});
% writetable(T,'D:\UQ_PCE\Diss\figures\counterexample_hinf_sample16.csv')
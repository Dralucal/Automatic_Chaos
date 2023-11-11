% Code related to Automatica Paper "On the Application of Galerkin Projection based Polynomial Chaos in Linear Systems and Control"
% by LL Evangelisti and H Pfifer
clear all; clc; close all;

%% 
load('projPCECoEx20.mat') % this matrix contains the projected DeltaPi matrix computed by the PolyChaos.jl Julia Toolbox

p = ureal('p',0,'Range',[-1,1]);
A = 0.01*[128*p^2-72*p-32,  295*p^2-199*p+4,   165*p^2-234*p+46; ...
		  -82*p^2-59*p+270, -266*p^2+144*p-73,  -147*p^2-210*p+286; ...
		  70*p^2+296*p-80,  43*p^2+96*p+8,      15*p^2+146*p-251];
% robstab(ss(A,[],[],[]))
[Me,Deltae,BLKSTRUCT] = lftdata(A);
Deltae
nd = length(Deltae);
M11 = Me(1:nd,1:nd); M12 = Me(1:nd,(nd+1):end); M21 = Me((nd+1):end,1:nd); M22 = Me((nd+1):end,(nd+1):end);
% save('D:\UQ_PCE\PCE\MA.mat','M11','M12','M21','M22','-v7')
%%
figure(1); 
maxdeg=16;
degvec = 0:1:maxdeg;
cm = jet(length(degvec)); set(groot,'defaultAxesColorOrder',cm);
x0 = [-0.7532899339782757; ...
      -0.6522776674778026;...
      0.08419097265575552];
nt=100;
t = linspace(0,15,nt)';
for deg = degvec
    L = deg + 1;
    e1 = zeros(L,1); e1(1) = 1;
    DeltaPiRep = kron(eye(BLKSTRUCT(1).Occurrences), DeltaPi(1:L,1:L));
    Api = lft( DeltaPiRep, kron(Me, eye(L)) );
    Xmean = zeros(nt,3);
    for i = 1:nt
       Xi = expm(Api*t(i))*kron(x0,e1); 
       Xmean(i,:) = Xi(1:L:(3*L));
    end
    figure(1);
    subplot(3,1,1)
    plot(t,Xmean(:,1)); ylabel('x_1'); xlabel('t')
    hold on
    [deltavec, weightvec] = legzo(L);
    Gauss_mean = zeros(nt,1);
    Gauss_realizations = zeros(nt,L);
    for i = 1:L
        As = usubs(A,'p',deltavec(i));
        for j = 1:nt
            xj = expm(As*t(j))*x0;
            Gauss_realizations(j,i) = xj(1);
        end
    end
    Gauss_mean = (Gauss_realizations * weightvec')./2;
    figure(2); plot(t, Gauss_mean - Xmean(:,1) ); title('Error between Gauss Quadrature and PCE')
    max(abs(Gauss_mean - Xmean(:,1)))
    %ha = gca;
    %ha.ColorOrderIndex = deg+1;
    figure(1); subplot(3,1,2)
    plot(t,Xmean(:,2)); ylabel('x_2'); xlabel('t')
    hold on
    %ha.ColorOrderIndex = deg+1;
    subplot(3,1,3)
    plot(t,Xmean(:,3)); ylabel('x_3'); xlabel('t')
    hold on
end
%%
% maxdeg=16;
% degvec = 0:1:maxdeg;
% cm = jet(length(degvec)); set(groot,'defaultAxesColorOrder',cm);
% x0 = [-0.7532899339782757; ...
%       -0.6522776674778026;...
%       0.08419097265575552];
% nt=100;
% t = linspace(0,25,nt)';
% X1mean = zeros(nt,maxdeg+1);
% for deg = degvec
%     L = deg + 1;
%     e1 = zeros(L,1); e1(1) = 1;
%     DeltaPiRep = kron(eye(BLKSTRUCT(1).Occurrences), DeltaPi(1:L,1:L));
%     Api = lft( DeltaPiRep, kron(Me, eye(L)) );
%     for i = 1:nt
%        Xi = expm(Api*t(i))*kron(x0,e1); 
%        Xmean(i,:) = Xi(1:L:(3*L));
%     end
%     X1mean(:,deg+1) = Xmean(:,1);
% end
% %%
% T = table(t,X1mean,'VariableNames',{'t','X1mean'});
% writetable(T,'D:\UQ_PCE\CDC23\figures\counterexample_GQ_sim.csv')
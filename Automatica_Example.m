% Code for the Plots to Automatica Paper "On the Application of Galerkin Projection based Polynomial Chaos in Linear Systems and Control"
% by LL Evangelisti and H Pfifer
clear all; clc;

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
N=100;
As = usample(A,N);
eigsmpl = NaN(3,N);
for i = 1:N
   eigsmpl(:,i) = eig(As(:,:,i));
end
%% N=100 randomly picked realizations are in blue
eigsmpl = reshape(eigsmpl,1,N*3);
figure; plot(real(eigsmpl),imag(eigsmpl),'bx'); hold on
% k = convhull(real(eigsmpl)',imag(eigsmpl)');
% plot(real(eigsmpl(k)),imag(eigsmpl(k)),'r-')
%% eigenvalues of Galerkin projected LFT in black are realizations
load('CALCdeltaPi10.mat') % this matrix contains the projected DeltaPi matrix computed by the PolyChaos.jl Julia Toolbox

maxdeg=10;
for deg = 0:maxdeg
    L = deg + 1;
    DeltaPiRep = kron(eye(BLKSTRUCT(1).Occurrences), DeltaPi(1:L,1:L));
    Api = lft( DeltaPiRep, kron(Me, eye(L)) );
    eigpi = eig(Api);
    plot(real(eigpi),imag(eigpi),'kx');
%     any(real(eig(Api)) >= 0);
%     e1 = zeros(L,1); e1(1) = 1;
% %     hinfnorm(ss(Api,kron(eye(3),e1),kron(eye(3),e1'),[]))
% 
%     gTry = 31; gLow = 0; gUpp = Inf; RelTol = 1e-2; AbsTol = 1e-3;
% 
%     while( (gUpp - gLow >= RelTol*gUpp) && (gUpp - gLow >= AbsTol) )
%         Ham_pi = [Api, eye(L*3)/gTry^2; -eye(L*3) -Api'];
%         if(min(abs(real(eig(Ham_pi))))<1e-3)
%     %     if(sol.x(1) < T)
%             gLow = gTry;
%             if(isfinite(gUpp))
%                 gTry = (gLow + gUpp)/2;
%             else
%                 gTry = gTry*2;
%             end
%         else
%             gUpp = gTry;
%             gTry = (gLow + gUpp)/2;
%         end
%     end
%     
%     gUpp
end
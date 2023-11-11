% Code for the Plots to Automatica Paper "On the Application of Galerkin Projection based Polynomial Chaos in Linear Systems and Control"
% by LL Evangelisti and H Pfifer
clear all; clc; close all;

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
N=5e1;
%As = usample(A,N);
deltavec = linspace(-1,1,N);
As=squeeze( usubs(A,'p',deltavec) );
eigsmpl = NaN(3,N);
figure; cm = jet(N);
set(groot,'defaultAxesColorOrder',cm);
for i = 1:N
   eigsmpl(:,i) = eig(As(:,:,i));
   plot(real(eigsmpl(:,i)),imag(eigsmpl(:,i)),'x');
   hold on;
end
%%
% plot(real(eigsmpl(1,:)),imag(eigsmpl(1,:)))
% hold on;
% plot(real(eigsmpl(2,:)),imag(eigsmpl(2,:)))
% plot(real(eigsmpl(3,:)),imag(eigsmpl(3,:)))
%%
% T = table(real(eigsmpl(1,:))',imag(eigsmpl(1,:))',...
%           real(eigsmpl(2,:))',imag(eigsmpl(2,:))',...
%           real(eigsmpl(3,:))',imag(eigsmpl(3,:))',...
%           deltavec',...
%           'VariableNames',{'real1','imag1','real2','imag2','real3','imag3','delta'});
% writetable(T,'D:\UQ_PCE\CDC23\figures\counterexample_heatmap_data.csv')
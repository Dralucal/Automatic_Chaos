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
N=1e3;
%As = usample(A,N);
deltavec = linspace(-1,1,N);
As=squeeze( usubs(A,'p',deltavec) );
eigsmpl = NaN(3,N);
figure; cm = jet(N);
set(groot,'defaultAxesColorOrder',cm);
for i = 1:N
   eigsmpl(:,i) = sort(eig(As(:,:,i)),'ComparisonMethod','real');
   plot(real(eigsmpl(:,i)),imag(eigsmpl(:,i)),'x');
   hold on;
end
%%
% plot(real(eigsmpl(1,:)),imag(eigsmpl(1,:)))
% hold on;
% plot(real(eigsmpl(2,:)),imag(eigsmpl(2,:)))
% plot(real(eigsmpl(3,:)),imag(eigsmpl(3,:)))
%%
ireal = abs(imag( (eigsmpl(2,:)) )) < 1e-2;
ifirstreal=find(ireal,1,'first');
ifirstbranch = 1:(ifirstreal-1);
ilastreal=find(ireal,1,'last');
isecbranch = ifirstreal : ilastreal;
ithirdbranch = (ilastreal+1):N;
% T = table(real(eigsmpl(1,:))',imag(eigsmpl(1,:))',deltavec',...
%           'VariableNames',{'real1','imag1','delta1'});
% writetable(T,'D:\UQ_PCE\CDC23\figures\counterexample_eigvalue1_comparison.csv')
% Error using table (line 231)
% % All input variables must have the same number of rows.
% T = table(real(eigsmpl(2,ifirstbranch))',imag(eigsmpl(2,ifirstbranch))',deltavec(ifirstbranch)',...
%           real(eigsmpl(3,ifirstbranch))',imag(eigsmpl(3,ifirstbranch))',deltavec(ifirstbranch)',...
%           'VariableNames',{ 'real21','imag21','delta21',...
%                             'real31','imag31','delta31'});
% writetable(T,'D:\UQ_PCE\CDC23\figures\counterexample_eigvalue231_comparison.csv')
% 
% T = table(real(eigsmpl(2,isecbranch))',imag(eigsmpl(2,isecbranch))',deltavec(isecbranch)',...
%           real(eigsmpl(3,isecbranch))',imag(eigsmpl(3,isecbranch))',deltavec(isecbranch)',...
%           'VariableNames',{ 'real22','imag22','delta22',...
%                             'real32','imag32','delta32'});
% writetable(T,'D:\UQ_PCE\CDC23\figures\counterexample_eigvalue232_comparison.csv')
% 
% T = table(real(eigsmpl(2,ithirdbranch))',imag(eigsmpl(2,ithirdbranch))',deltavec(ithirdbranch)',...
%           real(eigsmpl(3,ithirdbranch))',imag(eigsmpl(3,ithirdbranch))',deltavec(ithirdbranch)',...
%           'VariableNames',{ 'real23','imag23','delta23',...
%                             'real33','imag33','delta33'});
% writetable(T,'D:\UQ_PCE\CDC23\figures\counterexample_eigvalue233_comparison.csv')
%%
load('projPCECoEx20.mat') % this matrix contains the projected DeltaPi matrix computed by the PolyChaos.jl Julia Toolbox
figure

deg = 0; L = deg + 1;
DeltaPiRep = kron(eye(BLKSTRUCT(1).Occurrences), DeltaPi(1:L,1:L));
Api = lft( DeltaPiRep, kron(Me, eye(L)) );
eigpi0 = eig(Api);
plot(real(eigpi0),imag(eigpi0),'rsquare');

hold on

deg = 8; L = deg + 1;
DeltaPiRep = kron(eye(BLKSTRUCT(1).Occurrences), DeltaPi(1:L,1:L));
Api = lft( DeltaPiRep, kron(Me, eye(L)) );
eigpi8 = eig(Api);
plot(real(eigpi8),imag(eigpi8),'k*');

deg = 16; L = deg + 1;
DeltaPiRep = kron(eye(BLKSTRUCT(1).Occurrences), DeltaPi(1:L,1:L));
Api = lft( DeltaPiRep, kron(Me, eye(L)) );
eigpi16 = eig(Api);
plot(real(eigpi16),imag(eigpi16),'kx');
%%
% T = table(real(eigpi0),imag(eigpi0),'VariableNames',{ 'real','imag'});
% writetable(T,'D:\UQ_PCE\CDC23\figures\counterexample_eigvalue0_LFT.csv')
% T = table(real(eigpi8),imag(eigpi8),'VariableNames',{ 'real','imag'});
% writetable(T,'D:\UQ_PCE\CDC23\figures\counterexample_eigvalue8_LFT.csv')
% T = table(real(eigpi16),imag(eigpi16),'VariableNames',{ 'real','imag'});
% writetable(T,'D:\UQ_PCE\CDC23\figures\counterexample_eigvalue16_LFT.csv');
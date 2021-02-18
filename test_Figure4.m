clear all, clc, close all
load('DT2T120201208.mat')

%%
T1list = logspace(-1,1,20);
T2list = logspace(-1,1,20);
DifCoe = linspace(0,10,20);
D2D = squeeze( D3D(:,:,1) );
% y = D2D(:);
y =  D3D(:);
y = y(:);
figure,plot(y);

K3 = -(1-2*exp(-taut1.'*(1./T1list)));
K2 = exp(-taut2.'*(1./T2list));
K1 = exp(-b.'*DifCoe);
% ------ SVD Compress --------------------
K = kron(K3,kron(K2,K1)) ;
[u,sv,v]=svd(K,'econ');
figure,plot(diag(sv),'ro');
u = u(:,1:50);
v = v(:,1:50);
s = u'*y;
K = u'*K; 
% ---------------------------------------------
lambda =  0.1;
    rel_tol = 1e-8;
    sn = s; 
    sn_max = max(sn);
    sn = sn./sn_max;
    tic
    [x,status]=l1_ls_nonneg( K,sn,lambda,rel_tol);
    toc
    x =x./max(x); 
    X_tnipm = reshape( x,size(K1,2),size(K2,2),[]);
% figure,plot(x);   
%% Show Results
% yfit = K*x;
% figure,plot(T1list,x);
% figure,hold on,plot(taut1,y,'ro');
% plot(taut1,yfit)

% figure,
% contour(T2list,T1list,X_tnipm),set(gca,'XScale','log','YScale','log');

%
Dlable = DifCoe;%
T2lable = T2list;
T1lable = T1list;

X3D_tnipm = X_tnipm;
X3D_tnipm(X3D_tnipm < 0.2) = 0; % thresholding for the clear result show
figure,hold on,contourslice(Dlable,T2lable,T1lable,X3D_tnipm,Dlable,T2lable,T1lable,10); %DifCoe,T2list,T1list,
xlabel('D');ylabel('T2');zlabel('T1');
set (gca,'XGrid','on','YGrid','on','ZGrid','on',  'YTick',[0.1,1,10],'ZTick',[0.1,1,10])
set(gca,'YScale','log','ZScale','log')
% title('3D Laplace Spectra');
view(3); axis tight
xlim([min(DifCoe),max(DifCoe)]),ylim([min(T2list),max(T2list)]),zlim([min(T1list),max(T1list)])

close all
clear all
load('T2T2_50ms.mat')


   %  --------------------- TNIPM T2-T2 2D Reconstruction --------------------------------
    lambda =  0.0001;       %
    rel_tol = 1e-4;
    
    T21 = logspace(log10(1),log10(100),32);
    T22 = logspace(log10(1),log10(100),32);
    KT21 = exp(1000*t1values*(-1./T21));
    KT22 = exp(1000*t2values*(-1./T22));
    K_KT2T2 = kron(KT22,KT21);
    
    sn =  TT2data; 
    sn = sn(:);
    [u,s,v] = svd(K_KT2T2,'econ');
    u = u(:,1:50);
    v = v(:,1:50);
    sn = u'*sn;
    K = u'*K_KT2T2; 
    
    
    sn_max = max(sn);
    sn = sn./sn_max;
    [x_T2T2,status]=l1_ls_nonneg( K,sn,lambda,rel_tol);
    x_T2T2 = sn_max*x_T2T2;
    X_tnipm_T2T2 = reshape( x_T2T2,size(KT21,2),size(KT22,2) );
    
    figure,contour(T21,T22,X_tnipm_T2T2,200),colorbar
    title('T2-T2 via TNIPM');
    xlabel('T2£¨s£©'),ylabel('T2£¨s£©');
    set(gca,'XScale','log','YScale','log');
    set(gca,'YTick', [],'XTick', []);
    xlim([1.5,90]),ylim([1.5,90])
    figure,mesh(T21,T22,X_tnipm_T2T2)
     set(gca,'XScale','log','YScale','log');
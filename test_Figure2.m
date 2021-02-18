clear all
load('HP.mat')
%%  ----------------------- Initialization --------------------------------

rel_tol = 1e-5;
T2 = T2(2:2:90);
KT2 = exp( -t2list*(1./T2).' );
for iter=1:size(KT2,2)
   KT2_fa(iter) = norm(KT2(:,iter),2); 
end    
KT2_normal = KT2*diag(1./KT2_fa);

    rel_tol = 1e-8;
    D = linspace(1,20,45)*1e-9;
    KD = exp( -(qsqu).'*D );
    for iter=1:size(KD,2)
       KD_fa(iter) = norm(KD(:,iter),2); 
    end    
    KD_normal = KD*diag(1./KD_fa);

   %%  --------------------- TNIPM D-T2 2D Reconstruction without SVD Compress --------------------------------
    lambda =  1;
    rel_tol = 1e-8;
    K_DT2 = kron(KT2,KD);
    sn =  Dt2scal.' ; 
    sn = sn(:);
    sn_max = max(sn);
    sn = sn./sn_max;
    [x_DT2,status]=l1_ls_nonneg( K_DT2,sn,lambda,rel_tol);
    x_DT2 =x_DT2./max(x_DT2); 
    X_tnipm_DT2 = reshape( x_DT2,size(KD,2),size(KT2,2) );
    
    figure,contour(T2,10^8.4*D,X_tnipm_DT2,30),colorbar
    title('HP DT2 Spectra via TNIPM');
    xlabel('T2£¨s£©'),ylabel('D(10^-^9m^2s^-^10');
    set(gca,'XScale','log','YScale','log');
   xlim([0.1,10]),ylim([0.1,10]),
    set(gca,'YTick', [],'XTick', []);

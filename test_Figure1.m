clear all,close all,clc 
load('M6.mat') ;  %%% 
 

 
%---------------Initialization-------------------- 
dt=real(NmrData.SPECTRA); % processed experimental data
  if size(dt,1)>size(dt,2)
     dt = dt.';  
  end  
g=100*NmrData.Gzlvl; % gradient values
BD=NmrData.DELTAOriginal; % diffusion time
LD=NmrData.deltaOriginal; % diffusion encoding time
cs=NmrData.Specscale;     % chemical shift
gamma=4257.7; 
g2=(2*pi*gamma*g*LD).^2*(BD-LD/3)*1e4;
g2 = g2*1e-10;

cst2=[0.5 5.5];
% the chemical shift range of the reference peak
for i=1:size(cs,2)
 if(cs(i)>cst2(2))  % default: the value of chemical shift increases with i.
    cst=i;
    break;
 end
end
for i=1:size(cs,2)
 if(cs(i)>cst2(1))
    cfh=i;
    break;
 end
end
cdt = dt(:,cfh:cst);
cdt = cdt./max(max(cdt));
ccs = cs(cfh:cst);
Spec = cdt(1,:);
figure,plot(Spec);set(gca,'Xdir','reverse');%xlim([cst2]);

PeaBou = [574,690,1042,1176,1295,1429,1481,1597,1625,...
    1702,1950,2032,2225,2351,4052,4169,4607,4668,4792,...
    4883,4965,5015,5033,5097,5134,5195,5294,5346,5397,...
    5487,5514,5556,5768,5828,6046,6101,6109,6176,7989,8045];
PeaNum = length(PeaBou)/2;
figure,hold on
for it = 1:PeaNum
    PeaDat(:,it) = sum(cdt(:,PeaBou(it*2-1):PeaBou(it*2)),2);
    plot(PeaDat(:,it));
end



aa = 0:0.1:20;
K = exp(-g2.'*aa);

for iter=1:size(K,2)
   k_fa(iter) = norm(K(:,iter),2); 
end    
K_normal = K*diag(1./k_fa);
lambda =[0.1,0.1,0.05,0.05,0.005,0.1,0.1,0.1,0.1,0.1,0.1,0.01,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1];% individual regularization parameter setting
rel_tol = 1e-8;
X_Dosy = zeros(length(aa),length(Spec));
%--------------------- TNIPM --------------------------------
for it = 1:PeaNum
    
    sn = PeaDat(:,it);
    sn_max = max(sn);
    sn = sn./sn_max;
     
   [x_tnipm,status]=l1_ls_nonneg(K_normal,sn,lambda(it),rel_tol);
    x_tnipm = diag(1./k_fa)*x_tnipm*sn_max;
     x_tnipm =x_tnipm./max(x_tnipm);
    X_Dosy(:,PeaBou(2*it-1):PeaBou(2*it)) = repmat(x_tnipm,1,PeaBou(2*it)-PeaBou(2*it-1)+1)*diag(Spec(PeaBou(2*it-1):PeaBou(2*it)));
end
%%
figure,
      ax1 = axes('position',[0.1 0.7 0.8 0.3]);
      plot(ax1,ccs,Spec);set(gca,'Xdir','reverse');axis off;
      xlim([0.5,5.5]);
      ax2 = axes('position',[0.1 0.23 0.8 0.5]);
      contour(ax2,ccs,aa,X_Dosy,40);xlabel('Chemical Shift/ppm');ylabel('Diffusion Coefficient/10^-^1^0m^2s^-^1');
      set(ax2,'Ydir','reverse','Xdir','reverse'); 
  
      Comp_pos = [3.2,4.3,5.0,6.3,8.2,10.6];
      Comp_num = length(Comp_pos);
      set(ax2,'YTick', Comp_pos);
      xlim([0.5,5.5]);ylim([2,13]);
      set(gca,'FontWeight','bold','FontSize',14);
      for i = 1:Comp_num
           line(get(ax2,'xlim'),[Comp_pos(i),Comp_pos(i)],'LineWidth',0.8,'color',[0.75 0.75 0.75],'LineStyle','--' );
      end  
 
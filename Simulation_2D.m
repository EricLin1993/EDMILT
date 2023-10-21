clear all, clc, close all

t1 = linspace(1,10,10) ;
t2 = linspace(1,10,10) ;
 
N1 = length(t1);
N2 = length(t2);


alpha1 = -[0.2,0.1,0.6, 0.15,0.05,0.7,0.05,0.25,0.85,0.01,0.01,0.85,0.25];
alpha2 = -[0.2,0.6,0.1,0.05,0.05,0.05,0.7,0.85,0.25,0.85,0.25,0.01,0.01];


Y = 0;
for it = 1:4  % Select 4 Components for Simulation
   Y = Y + exp(alpha1(it)*t1.')*exp(alpha2(it)*t2);
end

%%  
a = linspace(0,1.2,25) ;
l1 = length(a);
l2 = length(a);
l3 = length(a);

K1 = exp(-t1.'*a);
K2 = exp(-t2.'*a);

% ----------  2D simulation  ----------------
lambda =  0.01;
 K{1}=K1;   
 K{2}=K2;  
    tic
    [ X ] = EDMILT( K,Y,lambda);
    toc
figure,contour(X);






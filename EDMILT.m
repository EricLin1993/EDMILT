function [ X ] = EDMILT( K,Y,lambda)
%% Author: Enping Lin
%   Date : 2021.1
%    
%  Input: K, cell array, multi-Laplace kernel
%         Y, multidimensional acquired data
%         lambda, regularized parameter for l1 norm
%  Output: X : reconstructed multidimensional spectra
%

     tic;
     kn = length(K);
     if kn ~= ndims(Y)
         error('K,Y did not agree in dimension') ;
     end  
     A = 1;
     ln =1;
     for it =1:kn  
         A = kron(K{it},A);
         xdn(it) = size(K{it},2);
         ln = ln*size(K{it},1);
     end
     y = Y(:);
     if length(y) ~= ln
         error('K, Y did not agree in length');
     end    
     if kn ~= ndims(Y)
         error('K,Y did not agree in dimension') 
     end  
     % ---------- SVD Compressed ----------------------
     [u,sv,v]=svd(A,'econ');
     svv = diag(sv);
     svv =svv/max(svv);
     svdn = sum(svv>0.001);
     figure,plot(diag(sv),'ro');
     u = u(:,1:svdn);
     y = u'*y;
     A = u'*A; 
     %-----------------Iterative Optimized Algorithm --------------------------
     % ------------------------------------------------------------------------
     y = y./max(y);
     [z] = TNIPM_ILT(A,y,lambda); % from Kim et al.'s code with some simplification for Multidimensional ILT
     X = reshape(z,xdn);
     EDMILT_TIME = toc
end


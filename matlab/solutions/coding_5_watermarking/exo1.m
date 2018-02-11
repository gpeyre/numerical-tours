rho = sqrt(P)/norm(w.*abs(x0)) * 10^( -psnr_embedding/20 );
disp(['rho = ' num2str(rho,3) '.']);

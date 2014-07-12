nb = length(y);
e1 = nb/n; % number of bit per symbol
% comparison with entropy bound
e = -sum( log2(h).*h ); % you have to use here the formula of the entropy
disp( strcat(['Entropy=' num2str(e, 3) ', arithmetic=' num2str(e1,3) '.']) );

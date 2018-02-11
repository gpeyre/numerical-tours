function y = load_signal(name, n, options)

% load_signal - load a 1D signal
%
%   y = load_signal(name, n, options);
%
% name is a string that can be :
%   'regular' (options.alpha gives regularity)
%   'step', 'rand',
%   'gaussiannoise' (options.sigma gives width of filtering in pixels),
%   [natural signals]
%   'tiger', 'bell', 'bird'
%   [WAVELAB signals]
%   'HeaviSine', 'Bumps', 'Blocks',
%   'Doppler', 'Ramp', 'Cusp', 'Sing', 'HiSine',
%   'LoSine', 'LinChirp', 'TwoChirp', 'QuadChirp',
%   'MishMash', 'WernerSorrows' (Heisenberg),
%   'Leopold' (Kronecker), 'Piece-Regular' (Piece-Wise Smooth),
%	'Riemann','HypChirps','LinChirps', 'Chirps', 'Gabor'
%	'sineoneoverx','Cusp2','SmoothCusp','Gaussian'
%	'Piece-Polynomial' (Piece-Wise 3rd degree polynomial)

if nargin<2
    n = 1024;
end
options.null = 0;
if isfield(options, 'alpha')
    alpha = options.alpha;
else
    alpha = 2;
end

options.rep = '';
switch lower(name)
    case 'regular'
        y = gen_signal(n,alpha);
    case 'step'
        y = linspace(0,1,n)>0.5;
    case 'stepregular'
        y = linspace(0,1,n)>0.5; y=y(:);
        a = gen_signal(n,2); a = a(:);
        a = rescale(a,-0.1,0.1);
        y = y+a;
    case 'gaussiannoise'
        % filtered gaussian noise
        y = randn(n,1);
        if isfield(options, 'sigma')
            sigma = options.sigma; % variance in number of pixels
        else
            sigma = 20;
        end
        m = min(n, 6*round(sigma/2)+1);
        h = compute_gaussian_filter(m,sigma/(4*n),n);
        options.bound = 'per';
        y = perform_convolution(y,h, options);
    case 'rand'
        if isfield(options, 'p1')
            p1 = options.p1;
        else
            c = 10;
            p1 = 1:c; p1 = p1/sum(p1);
        end
        p1 = p1(:); c = length(p1);
        if isfield(options, 'p2')
            p2 = options.p2;
        else
            if isfield(options, 'evol')
                evol = options.evol;
            else
                evol = 0;
            end
            p2 = p1(:) + evol*(rand(c,1)-0.5);
            p2 = max(p2,0); p2 = p2/sum(p2);
        end
        y = zeros(n,1);
        for i=1:n
            a = (i-1)/(n-1);
            p = a*p1+(1-a)*p2; p = p/sum(p);
            y(i) = rand_discr(p, 1);
        end
    case 'bird'
        [y,fs] = load_sound([name '.wav'], n, options);
    case 'tiger'
        [y,fs] = load_sound([name '.au'], n, options);
    case 'bell'
        [y,fs] = load_sound([name '.wav'], n, options);
    otherwise
        y = MakeSignal(name,n);
end

y = y(:);

function y = gen_signal(n,alpha)

%   gen_signal -  generate a 1D C^\alpha signal of length n.
%
%   y = gen_signal(n,alpha);
%
%   The signal is scaled in [0,1].
%   
%   Copyright (c) 2003 Gabriel Peyr?

if nargin<2
    alpha = 2;
end

y = randn(n,1); 
fy = fft(y);
fy = fftshift(fy);
% filter with |omega|^{-\alpha}
h = (-n/2+1):(n/2);
h = (abs(h)+1).^(-alpha-0.5);
fy = fy.*h';
fy = fftshift(fy);
y = real( ifft(fy) );

y = (y-min(y))/(max(y)-min(y));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function sig = MakeSignal(Name,n)
% MakeSignal -- Make artificial signal
%  Usage
%    sig = MakeSignal(Name,n)
%  Inputs
%    Name   string: 'HeaviSine', 'Bumps', 'Blocks',
%            'Doppler', 'Ramp', 'Cusp', 'Sing', 'HiSine',
%            'LoSine', 'LinChirp', 'TwoChirp', 'QuadChirp',
%            'MishMash', 'WernerSorrows' (Heisenberg),
%            'Leopold' (Kronecker), 'Piece-Regular' (Piece-Wise Smooth),
%	     'Riemann','HypChirps','LinChirps', 'Chirps', 'Gabor'
%	     'sineoneoverx','Cusp2','SmoothCusp','Gaussian'
%	     'Piece-Polynomial' (Piece-Wise 3rd degree polynomial)
%    n      desired signal length
%  Outputs
%    sig    1-d signal
%
%  References
%    Various articles of D.L. Donoho and I.M. Johnstone
%
if nargin > 1,
    t = (1:n) ./n;
end
Name = lower(Name);
if strcmp(Name,'heavisine'),
    sig = 4.*sin(4*pi.*t);
    sig = sig - sign(t - .3) - sign(.72 - t);
elseif strcmp(Name,'bumps'),
    pos = [ .1 .13 .15 .23 .25 .40 .44 .65  .76 .78 .81];
    hgt = [ 4  5   3   4  5  4.2 2.1 4.3  3.1 5.1 4.2];
    wth = [.005 .005 .006 .01 .01 .03 .01 .01  .005 .008 .005];
    sig = zeros(size(t));
    for j =1:length(pos)
        sig = sig + hgt(j)./( 1 + abs((t - pos(j))./wth(j))).^4;
    end
elseif strcmp(Name,'blocks'),
    pos = [ .1 .13 .15 .23 .25 .40 .44 .65  .76 .78 .81];
    hgt = [4 (-5) 3 (-4) 5 (-4.2) 2.1 4.3  (-3.1) 2.1 (-4.2)];
    sig = zeros(size(t));
    for j=1:length(pos)
        sig = sig + (1 + sign(t-pos(j))).*(hgt(j)/2) ;
    end
elseif strcmp(Name,'doppler'),
    sig = sqrt(t.*(1-t)).*sin((2*pi*1.05) ./(t+.05));
elseif strcmp(Name,'ramp'),
    sig = t - (t >= .37);
elseif strcmp(Name,'cusp'),
    sig = sqrt(abs(t - .37));
elseif strcmp(Name,'sing'),
    k = floor(n * .37);
    sig = 1 ./abs(t - (k+.5)/n);
elseif strcmp(Name,'hisine'),
    sig = sin( pi * (n * .6902) .* t);
elseif strcmp(Name,'losine'),
    sig = sin( pi * (n * .3333) .* t);
elseif strcmp(Name,'linchirp'),
    sig = sin(pi .* t .* ((n .* .500) .* t));
elseif strcmp(Name,'twochirp'),
    sig = sin(pi .* t .* (n .* t)) + sin((pi/3) .* t .* (n .* t));
elseif strcmp(Name,'quadchirp'),
    sig = sin( (pi/3) .* t .* (n .* t.^2));
elseif strcmp(Name,'mishmash'),  % QuadChirp + LinChirp + HiSine
    sig = sin( (pi/3) .* t .* (n .* t.^2)) ;
    sig = sig +  sin( pi * (n * .6902) .* t);
    sig = sig +  sin(pi .* t .* (n .* .125 .* t));
elseif strcmp(Name,'wernersorrows'),
    sig = sin( pi .* t .* (n/2 .* t.^2)) ;
    sig = sig +  sin( pi * (n * .6902) .* t);
    sig = sig +  sin(pi .* t .* (n .* t));
    pos = [ .1 .13 .15 .23 .25 .40 .44 .65  .76 .78 .81];
    hgt = [ 4  5   3   4  5  4.2 2.1 4.3  3.1 5.1 4.2];
    wth = [.005 .005 .006 .01 .01 .03 .01 .01  .005 .008 .005];
    for j =1:length(pos)
        sig = sig + hgt(j)./( 1 + abs((t - pos(j))./wth(j))).^4;
    end
elseif strcmp(Name,'leopold'),
    sig = (t == floor(.37 * n)/n);  % Kronecker
elseif strcmp(Name,'riemann'),
    sqn = round(sqrt(n));
    sig = t .* 0;  % Riemann's Non-differentiable Function
    sig((1:sqn).^2) = 1. ./ (1:sqn);
    sig = real(ifft(sig));
elseif strcmp(Name,'hypchirps'), % Hyperbolic Chirps of Mallat's book
    alpha	= 15*n*pi/1024;
    beta    = 5*n*pi/1024;
    t  	= (1.001:1:n+.001)./n;
    f1      = zeros(1,n);
    f2      = zeros(1,n);
    f1  	= sin(alpha./(.8-t)).*(0.1<t).*(t<0.68);
    f2  	= sin(beta./(.8-t)).*(0.1<t).*(t<0.75);
    M  	= round(0.65*n);
    P 	= floor(M/4);
    enveloppe = ones(1,M); % the rising cutoff function
    enveloppe(1:P) = (1+sin(-pi/2+((1:P)-ones(1,P))./(P-1)*pi))/2;
    enveloppe(M-P+1:M) = reverse(enveloppe(1:P));
    env 	= zeros(1,n);
    env(ceil(n/10):M+ceil(n/10)-1) = enveloppe(1:M);
    sig     = (f1+f2).*env;
elseif strcmp(Name,'linchirps'), % Linear Chirps of Mallat's book
    b 	= 100*n*pi/1024;
    a 	= 250*n*pi/1024;
    t 	= (1:n)./n;
    A1 	= sqrt((t-1/n).*(1-t));
    sig	= A1.*(cos((a*(t).^2)) + cos((b*t+a*(t).^2)));
elseif strcmp(Name,'chirps'), % Mixture of Chirps of Mallat's book
    t 	= (1:n)./n.*10.*pi;
    f1 	= cos(t.^2*n/1024);
    a 	= 30*n/1024;
    t 	= (1:n)./n.*pi;
    f2 	= cos(a.*(t.^3));
    f2 	= reverse(f2);
    ix 	= (-n:n)./n.*20;
    g 	= exp(-ix.^2*4*n/1024);
    i1 	= (n/2+1:n/2+n);
    i2 	= (n/8+1:n/8+n);
    j  	= (1:n)/n;
    f3 	= g(i1).*cos(50.*pi.*j*n/1024);
    f4 	= g(i2).*cos(350.*pi.*j*n/1024);
    sig 	= f1+f2+f3+f4;
    enveloppe = ones(1,n); % the rising cutoff function
    enveloppe(1:n/8) = (1+sin(-pi/2+((1:n/8)-ones(1,n/8))./(n/8-1)*pi))/2;
    enveloppe(7*n/8+1:n) = reverse(enveloppe(1:n/8));
    sig 	= sig.*enveloppe;
elseif strcmp(Name,'gabor'), % two modulated Gabor functions in
    % Mallat's book
    N = 512;
    t = (-N:N)*5/N;
    j = (1:N)./N;
    g = exp(-t.^2*20);
    i1 = (2*N/4+1:2*N/4+N);
    i2 = (N/4+1:N/4+N);
    sig1 = 3*g(i1).*exp(i*N/16.*pi.*j);
    sig2 = 3*g(i2).*exp(i*N/4.*pi.*j);
    sig = sig1+sig2;
elseif strcmp(Name,'sineoneoverx'), % sin(1/x) in Mallat's book
    N = 1024;
    a = (-N+1:N);
    a(N) = 1/100;
    a = a./(N-1);
    sig = sin(1.5./(i));
    sig = sig(513:1536);
elseif strcmp(Name,'cusp2'),
    N = 64;
    a = (1:N)./N;
    x = (1-sqrt(a)) + a/2 -.5;
    M = 8*N;
    sig = zeros(1,M);
    sig(M-1.5.*N+1:M-.5*N) = x;
    sig(M-2.5*N+2:M-1.5.*N+1) = reverse(x);
    sig(3*N+1:3*N + N) = .5*ones(1,N);
elseif strcmp(Name,'smoothcusp'),
    sig = MakeSignal('Cusp2');
    N = 64;
    M = 8*N;
    t = (1:M)/M;
    sigma = 0.01;
    g = exp(-.5.*(abs(t-.5)./sigma).^2)./sigma./sqrt(2*pi);
    g = fftshift(g);
    sig2 = iconv(g',sig)'/M;
elseif strcmp(Name,'piece-regular'),
    sig1=-15*MakeSignal('Bumps',n);
    t = (1:fix(n/12)) ./fix(n/12);
    sig2=-exp(4*t);
    t = (1:fix(n/7)) ./fix(n/7);
    sig5=exp(4*t)-exp(4);
    t = (1:fix(n/3)) ./fix(n/3);
    sigma=6/40;
    sig6=-70*exp(-((t-1/2).*(t-1/2))/(2*sigma^2));
    sig(1:fix(n/7))= sig6(1:fix(n/7));
    sig((fix(n/7)+1):fix(n/5))=0.5*sig6((fix(n/7)+1):fix(n/5));
    sig((fix(n/5)+1):fix(n/3))=sig6((fix(n/5)+1):fix(n/3));
    sig((fix(n/3)+1):fix(n/2))=sig1((fix(n/3)+1):fix(n/2));
    sig((fix(n/2)+1):(fix(n/2)+fix(n/12)))=sig2;
    sig((fix(n/2)+2*fix(n/12)):-1:(fix(n/2)+fix(n/12)+1))=sig2;
    sig(fix(n/2)+2*fix(n/12)+fix(n/20)+1:(fix(n/2)+2*fix(n/12)+3*fix(n/20)))=...
        -ones(1,fix(n/2)+2*fix(n/12)+3*fix(n/20)-fix(n/2)-2*fix(n/12)-fix(n/20))*25;
    k=fix(n/2)+2*fix(n/12)+3*fix(n/20);
    sig((k+1):(k+fix(n/7)))=sig5;
    diff=n-5*fix(n/5);
    sig(5*fix(n/5)+1:n)=sig(diff:-1:1);
    % zero-mean
    bias=sum(sig)/n;
    sig=bias-sig;
elseif strcmp(Name,'piece-polynomial'),
    t = (1:fix(n/5)) ./fix(n/5);
    sig1=20*(t.^3+t.^2+4);
    sig3=40*(2.*t.^3+t) + 100;
    sig2=10.*t.^3 + 45;
    sig4=16*t.^2+8.*t+16;
    sig5=20*(t+4);
    sig6(1:fix(n/10))=ones(1,fix(n/10));
    sig6=sig6*20;
    sig(1:fix(n/5))=sig1;
    sig(2*fix(n/5):-1:(fix(n/5)+1))=sig2;
    sig((2*fix(n/5)+1):3*fix(n/5))=sig3;
    sig((3*fix(n/5)+1):4*fix(n/5))=sig4;
    sig((4*fix(n/5)+1):5*fix(n/5))=sig5(fix(n/5):-1:1);
    diff=n-5*fix(n/5);
    sig(5*fix(n/5)+1:n)=sig(diff:-1:1);
    %sig((fix(n/20)+1):(fix(n/20)+fix(n/10)))=-ones(1,fix(n/10))*20;
    sig((fix(n/20)+1):(fix(n/20)+fix(n/10)))=ones(1,fix(n/10))*10;
    sig((n-fix(n/10)+1):(n+fix(n/20)-fix(n/10)))=ones(1,fix(n/20))*150;
    % zero-mean
    bias=sum(sig)/n;
    sig=sig-bias;
elseif strcmp(Name,'gaussian'),
    sig=GWN(n,beta);
    g=zeros(1,n);
    lim=alpha*n;
    mult=pi/(2*alpha*n);
    g(1:lim)=(cos(mult*(1:lim))).^2;
    g((n/2+1):n)=g((n/2):-1:1);
    g = rnshift(g,n/2);
    g=g/norm(g);
    sig=iconv(g,sig);
else
    disp(sprintf('MakeSignal: I don*t recognize <<%s>>',Name))
    disp('Allowable Names are:')
    disp('HeaviSine'),
    disp('Bumps'),
    disp('Blocks'),
    disp('Doppler'),
    disp('Ramp'),
    disp('Cusp'),
    disp('Crease'),
    disp('Sing'),
    disp('HiSine'),
    disp('LoSine'),
    disp('LinChirp'),
    disp('TwoChirp'),
    disp('QuadChirp'),
    disp('MishMash'),
    disp('WernerSorrows'),
    disp('Leopold'),
    disp('Sing'),
    disp('HiSine'),
    disp('LoSine'),
    disp('LinChirp'),
    disp('TwoChirp'),
    disp('QuadChirp'),
    disp('MishMash'),
    disp('WernerSorrows'),
    disp('Leopold'),
    disp('Riemann'),
    disp('HypChirps'),
    disp('LinChirps'),
    disp('Chirps'),
    disp('sineoneoverx'),
    disp('Cusp2'),
    disp('SmoothCusp'),
    disp('Gabor'),
    disp('Piece-Regular');
    disp('Piece-Polynomial');
    disp('Gaussian');
end

%
% Originally made by David L. Donoho.
% Function has been enhanced.
    
    
%   
% Part of WaveLab Version 802
% Built Sunday, October 3, 1999 8:52:27 AM
% This is Copyrighted Material
% For Copying permissions see COPYING.m
% Comments? e-mail wavelab@stat.stanford.edu
%   
    


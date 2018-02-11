function y = load_signal(Name, n, options)

// load_signal - load a signal by name
//
//   y = load_signal(Name, n, options);
//
//	Copyright (c) 2008 Gabriel Peyre

t = (1:n) ./n;
Name = lower(Name);
pi = 3.1416;
if strcmp(Name,'heavisine'),
    sig = 4.*sin(4*pi.*t);
    sig = sig - sign(t - .3) - sign(.72 - t);
elseif strcmp(Name,'bumps'),
    pos = [ .1 .13 .15 .23 .25 .40 .44 .65  .76 .78 .81];
    hgt = [ 4  5   3   4  5  4.2 2.1 4.3  3.1 5.1 4.2];
    wth = [.005 .005 .006 .01 .01 .03 .01 .01  .005 .008 .005];
    sig = zeros(length(t),1); t = t(:);
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
elseif strcmp(Name,'mishmash'),  // QuadChirp + LinChirp + HiSine
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
    sig = (t == floor(.37 * n)/n);  // Kronecker
elseif strcmp(Name,'riemann'),
    sqn = round(sqrt(n));
    sig = t .* 0;  // Riemann's Non-differentiable Function
    sig((1:sqn).^2) = 1. ./ (1:sqn);
    sig = real(ifft(sig));
elseif strcmp(Name,'hypchirps'), // Hyperbolic Chirps of Mallat's book
    alpha	= 15*n*pi/1024;
    beta    = 5*n*pi/1024;
    t  	= (1.001:1:n+.001)./n;
    f1      = zeros(1,n);
    f2      = zeros(1,n);
    f1  	= sin(alpha./(.8-t)).*(0.1<t).*(t<0.68);
    f2  	= sin(beta./(.8-t)).*(0.1<t).*(t<0.75);
    M  	= round(0.65*n);
    P 	= floor(M/4);
    enveloppe = ones(1,M); // the rising cutoff function
    enveloppe(1:P) = (1+sin(-pi/2+((1:P)-ones(1,P))./(P-1)*pi))/2;
    enveloppe(M-P+1:M) = reverse(enveloppe(1:P));
    env 	= zeros(1,n);
    env(ceil(n/10):M+ceil(n/10)-1) = enveloppe(1:M);
    sig     = (f1+f2).*env;
elseif strcmp(Name,'linchirps'), // Linear Chirps of Mallat's book
    b 	= 100*n*pi/1024;
    a 	= 250*n*pi/1024;
    t 	= (1:n)./n;
    A1 	= sqrt((t-1/n).*(1-t));
    sig	= A1.*(cos((a*(t).^2)) + cos((b*t+a*(t).^2)));
elseif strcmp(Name,'chirps'), // Mixture of Chirps of Mallat's book
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
    enveloppe = ones(1,n); // the rising cutoff function
    enveloppe(1:n/8) = (1+sin(-pi/2+((1:n/8)-ones(1,n/8))./(n/8-1)*pi))/2;
    enveloppe(7*n/8+1:n) = reverse(enveloppe(1:n/8));
    sig 	= sig.*enveloppe;
elseif strcmp(Name,'gabor'), // two modulated Gabor functions in
    // Mallat's book
    N = 512;
    t = (-N:N)*5/N;
    j = (1:N)./N;
    g = exp(-t.^2*20);
    i1 = (2*N/4+1:2*N/4+N);
    i2 = (N/4+1:N/4+N);
    sig1 = 3*g(i1).*exp(i*N/16.*pi.*j);
    sig2 = 3*g(i2).*exp(i*N/4.*pi.*j);
    sig = sig1+sig2;
elseif strcmp(Name,'sineoneoverx'), // sin(1/x) in Mallat's book
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
    sig = load_signal('Cusp2');
    N = 64;
    M = 8*N;
    t = (1:M)/M;
    sigma = 0.01;
    g = exp(-.5.*(abs(t-.5)./sigma).^2)./sigma./sqrt(2*pi);
    g = fftshift(g);
    sig2 = iconv(g',sig)'/M;
elseif strcmp(Name,'piece-regular'),
    sig1=-15*load_signal('Bumps',n)';
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
    diff1=n-5*fix(n/5);
    sig(5*fix(n/5)+1:n)=sig(diff1:-1:1);
    // zero-mean
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
    diff1=n-5*fix(n/5);
    sig(5*fix(n/5)+1:n)=sig(diff1:-1:1);
    //sig((fix(n/20)+1):(fix(n/20)+fix(n/10)))=-ones(1,fix(n/10))*20;
    sig((fix(n/20)+1):(fix(n/20)+fix(n/10)))=ones(1,fix(n/10))*10;
    sig((n-fix(n/10)+1):(n+fix(n/20)-fix(n/10)))=ones(1,fix(n/20))*150;
    // zero-mean
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
	error('Unknown signal');
end

y = sig(:);

endfunction
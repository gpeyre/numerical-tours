%% Multiplicative Cascade Synthesis of Signals and Textures
% This numerical tour explores multifractal signal and texture synthesis.

%%
% This tour is written by Pierre Chainais and Gabriel Peyré.

%%
% The processes we deal with belong to the family of _Infinitely Divisible Cascades_ (IDC).
% Only the simulation of the subfamily of _Compound Poisson Cascades_ (CPC) is really simple
% to implement in 1D, 2D or even ND. Indeed, the synthesis of CPC can be understood as the 
% product of a random number of indicator function of balls (of a 1D segment or a 2D ball) 
% with randomized radius and randomized amplitude. 

%% 
% If the distribution of the amplitudes and radii is well chosen, this
% leads to the synthesis of a function that is the density of a
% positive *scale invariant* measure. More precisely, this measure is *multifractal*. 

%%
% To obtain the final measure signal/image, the simulated measure density is integrated
% or pseudo-integrated thanks to some scale invariant low-pass filtering in the Fourier domain.

%%
% The application of cascade for texture synthesis is detailed in 

%%
% _Infinitely divisible cascades to model the statistics of natural images_,
% P. Chainais,   IEEE Trans. on Pattern Analysis and Machine Intelligence,
% Vol. 29 no 12, Dec. 2007.

%%
% Visit the homepage of <http://www.isima.fr/~pchainai/PUB/software.html
% Pierre Chainais> for additional information and softwares.

perform_toolbox_installation('signal', 'general', 'graph');


%% 1D Multiplicative Cascades
% We consider here the synthesis of 1D signals using multiplicative
% compound Poisson cascades (CPC).

%%
% Duration of the simulation : [0, T]

T = 10;

%%
% Resolution of the cascade (the smallest scale at which details will be added).

rmin = 0.02;

%%
% Sampling period of the simulation. Number of points.

Delta_t = rmin/2;
n = T/Delta_t+1;

%%
% Number of multipliers Wi.

lambda = (1/rmin-1)*(T+1); % to ensure scale invariance, scales are distributed as 1/r^2.
N = round(lambda);
% N = poissrnd(N,1,1); % rigorously, N is a Poisson r.v. of expectation lambda

%%
% Time positions of the Poisson point process in the time/scale plane
% are uniformly distributed to ensure stationarity. Side effects are avoided by extending
% the cascade to [-1/2, 0] and [T T+1/2].

ti = -1/2+rand(1,N) * (T+1);
ti_1 = -1/2+rand(1,round(T+1))*(T+1);
ti = [ti_1 ti];   % for exact scale invariance

%%
% Scales ri of the points.
% Should be distributed according to dr/r^2, to have scale invariance.

umax = 1/rmin-1;
ui = [zeros(1,length(ti_1)) rand(1,N) * umax ]; % ui = 1/ri-1
ri = (1+ui).^(-1);


%%
% Display the Poisson point process in the scale-space plane.

figure(1)
clf
plot(ti, ri, '.')
axis([-1/2,T+1/2,0,1]);

%%
% Parameters for the law of multipliers Wi (Wi>0). 
% Here we choose a log-normal law. Another simple possible choice
% is to set Wi=2/3 for all the Wi.

sigma2 = 0.2;
mu = -sigma2/2;

%% 
% Condition of non-degeneracy:

if -(exp(2*(sigma2+mu))-1)<-1
   disp('Be careful ! This cascade will degenerate as rmin -> 0 !')
end

%%
% Random log-normal multipliers.
N = length(ti);   % the number of multipliers = number of time-scale points.
Wi = exp( randn(N,1)*sqrt(sigma2)+mu );

%%
% Positions along time axis.

t = linspace(0,T,n);

%% 
% Initialize the signal and normalize the measure.
H1 = 1 - exp(mu+sigma2/2);
f = ones(n,1) * exp(H1) / rmin^H1;

%%
% We show here the first step of the multiplicative cascade: 
% iterations are on the multipliers (ti,ri,Wi). 

%%
% Select the points in the cone of influence of |(ti(1),ri(1))|.

i = 1;
I = find(abs(t-ti(i))<=ri(i)/2); % t belongs to a disk centered on ti(i)

%%
% Perform the multiplication with the random multiplier.

f(I) = f(I) * Wi(i);

clf
plot(t,f)
axis([0 T 0 1.1*max(f)])



%EXO
%% Perform the cascade. Display intermediate steps.
f = ones(n,1) * exp(H1) / rmin^H1;
k = 0;
clf
for i=1:N
   I = find(abs(t-ti(i))<=ri(i)/2);
   f(I) = f(I) * Wi(i);    
   if i==10 || i==80 || i==round(.1*N) || i==round(.5*N)
       k = k+1;
       subplot(2,2,k);
       plot(t,f)
		axis([0 T 0 1.1*max(f)])
   end
end
%EXO

%% 
% Display the random measure.
figure(1)
clf
plot(t,f)
axis([0 T 0 1.1*max(f)])



%EXO
%% Compute several realization for the same log-normal parameters.
clf
for k=1:4
   Wi = exp( randn(N,1)*sqrt(sigma2)+mu );
   ti = -1/2+rand(1,N) * (T+1);
   ti_1 = -1/2+rand(1,round(T+1))*(T+1);
   ti = [ti_1 ti];   % for exact scale invariance
   umax = 1/rmin-1;
   ui = [zeros(1,length(ti_1)) rand(1,N) * umax ]; % ui = 1/ri-1
   ri = (1+ui).^(-1);
   f = ones(n,1)* exp(H1) / rmin^H1;
   for i=1:N
       I = find(abs(t-ti(i))<=ri(i)/2);
       f(I) = f(I) * Wi(i);
   end
   subplot(2,2,k);
   plot(t,f)
   axis([0 T 0 1.1*max(f)])
end
%EXO

%EXO
%% Compute realizations for different log-normal parameters |mu| and |sigma2|.
%% Use the same distribution of points.
sigma2list = [0.1 0.1 0.5 0.5];
mult = [1 10 1 10];
mulist = -mult.*sigma2list;
ti = -1/2+rand(1,N) * (T+1);
ti_1 = -1/2+rand(1,round(T+1))*(T+1);
ti = [ti_1 ti];   % for exact scale invariance
umax = 1/rmin-1;
ui = [zeros(1,length(ti_1)) rand(1,N) * umax ]; % ui = 1-1/ri^{1+beta}
ri = (1+ui).^(-1);
clf
for k=1:4
   mu = mulist(k);
   sigma2 = sigma2list(k);
   if -(exp(2*(sigma2+mu))-1)<-1  % Condition of non-degeneracy
       disp('Be careful ! This cascade will degenerate as rmin -> 0 !')
   end
   H1 = 1 - exp(mu+sigma2/2);
   Wi = exp( randn(N,1)*sqrt(sigma2)+mu );
   f = ones(n,1)* exp(H1) / rmin^H1;
   for i=1:N
       I = find(abs(t-ti(i))<=ri(i)/2);
       f(I) = f(I) * Wi(i);
   end
   subplot(2,2,k);
   plot(t,f)
   axis([0 T 0 1.1*max(f)]) 
   title(['mu=' num2str(mu) ', sigma^2=' num2str(sigma2)]);
end
%EXO

%% 2D Multiplicative Cascades
% To generate 2D cascade, one needs to throw points on a 3D scale space
% domain.

%% 
% Size of the image.

n = 128;

%%
% Minimum scale, should be roughly |1/n|.

rmin = 1/n;

%%
% Maximum scale.
Xmax = 1;
Ymax = 1;

%%
% Number of points in the cascade.

lambda = 2/pi*(1/rmin^2-1)*(Xmax+1)*(Ymax+1); % density g(r)dr=4/pi/r^3 dr
N = round(lambda); % should be a Poisson r.v. with expectation lambda.

%%
% Scale of the points.
% Should be distributed according to 1/height^3, to have scaling
% invariance.

umax = (1/rmin^2-1);          % u will be a uniform variable in [0 1/rmin^2-1]
ui = rand(1,N) * umax;
ri = (1/rmin^2-ui).^(-1/2);

%%
% Position of the points.

xi = -1/2 + rand(1,N) * (Xmax+1);
yi = -1/2 + rand(1,N) * (Ymax+1);

%%
% Display the points in the scale-space plane.

clf
h = plot3(xi, yi, ri, '.');
axis([-1/2,Xmax+1/2,-1/2,Ymax+1/2,0,1]);


%%
% Parameters for the log-normal law.

sigma2 = 0.08;
mu = -sigma2/2;

%%
% Random log-normal multipliers.

Wi = exp( randn(N,1)*sqrt(sigma2)+mu );

%%
% Position in the X/Y plane.
% We enlarge the square in order to be able to use periodic boundary
% conditons.

x = linspace(0,Xmax,n);
y = linspace(0,Ymax,n);
[X,Y]= meshgrid(x,y);

%% 
% Initialization and normalization of the image.
H1 = 1 - exp(mu+sigma2/2);
f = ones(n)/rmin^H1;

%%
% We give here the example of the first mutiplication.

%%
% Localization of the signal locations that are influenced by the 
% scale/space point indexed by |(xi(1),yi(1),ri(1))|.
% This corresponds to the intersection of the image plane and a cone of influence.

i = 1;
I = find( (X-xi(i)).^2+(Y-yi(i)).^2 <=ri(i)^2/4 );

%%
% Multiplication of the image with the random multiplier.

f(I) = f(I) * Wi(i); 



%EXO
%% Perform the full cascade, display intermediate steps.
f = ones(n)/rmin^H1;
k = 0;
clf
for i=1:N  
   %progressbar(i,N);
   I = find( (X-xi(i)).^2+(Y-yi(i)).^2 <=ri(i)^2/4 );
   f(I) = f(I) * Wi(i); 
   if i==30 || i==200 || i==round(N*.5) || i==round(N*.8)  
       k = k+1;
       subplot(2,2,k);
       imageplot(f);
   end
end
%EXO

%%
% Display the image. It corresponds to a 2D multi-fractal measure.
clf
imageplot(f)


%%
% To compute the final texture, we perform a *spectral integration*, which
% corresponds to a low pass filtering.

%%
% Fourier frequency localizations.

x = [0:n/2 -n/2+1:-1];
[U,V] = meshgrid(x,x); 
S = sqrt(U.^2+V.^2); 
S(1,1) = 1;

%%
% Exponent of integration.
alpha = .5;

%%
% Fourier domain integration.

F = real( ifft2( fft2(f)./S.^alpha ) );


%EXO
%% Compute the fractional integration for several values of alpha.
alphalist = [.1 .3 .6 1];
for i=1:4
   alpha = alphalist(i);
   F = real( ifft2( fft2(f)./S.^alpha ) );
   subplot(2,2,i);
   imageplot(F, ['alpha=' num2str(alpha)]);
end
%EXO

%EXO
%% Perform the  cascade for several log-normal parameters |mu| and |sigma2|.
sigma2list = [0.02 0.02 0.2 0.2];
mult = [1 10 1 1];
mulist = -mult.*sigma2list/2;
k = 0;
clf;
for k=1:4
   mu = mulist(k);
   sigma2 = sigma2list(k);
   if -(exp(2*(sigma2+mu))-1)<-1  % Condition of non-degeneracy
       disp('Be careful ! This cascade will degenerate as rmin -> 0 !')
   end
   H1 = 1 - exp(mu+sigma2/2);
   Wi = exp( randn(N,1)*sqrt(sigma2)+mu );
   f = ones(n)/rmin^H1;
   for i=1:N
       %progressbar(i,N);
       I = find( (X-xi(i)).^2+(Y-yi(i)).^2 <=ri(i)^2/4 );
       f(I) = f(I) * Wi(i);
   end
   subplot(2,2,k);
   imageplot(f, ['mu=' num2str(mu) ', sigma^2=' num2str(sigma2)]);
end
%EXO
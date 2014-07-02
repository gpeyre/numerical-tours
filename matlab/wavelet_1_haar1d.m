%% 1-D Haar Wavelets
% This numerical tour explores 1-D multiresolution analysis with Haar
% transform. It was introduced in 1910 by Haar  <#biblio [Haar1910]>
% and is arguably the first example of wavelet basis.

perform_toolbox_installation('signal', 'general');


%% Forward Haar Transform
% The Haar transform is the simplest orthogonal wavelet transform. It is
% computed by iterating difference and averaging between odd and even
% samples of the signal.

%%
% Size \(N\) of the signal.

N = 512;

%%
% First we load a 1-D signal \(f \in \RR^N\).

name = 'piece-regular';
f = rescale( load_signal(name, N) );

%%
% The Haar transform operates over \(J = \log_2(N)-1\) scales.
% It computes a series of coarse scale and fine scale coefficients
% \(a_j, d_j \in \RR^{N_j}\) where \(N_j=2^j\).

J = log2(N)-1;

%%
% The forward Haar transform computed \( \Hh(f) = (d_j)_{j=0}^J \cup \{a_0\} \).
% Note that the set of coarse scale coefficients \((a_j)_{j>0}\) are not
% stored.

%%
% This transform is orthogonal, meaning \( \Hh \circ \Hh^* = \text{Id} \),
% and that there is the following conservation of energy
% \[ \sum_i \abs{f_i}^2 = \norm{f}^2 = \norm{\Hh f}^2 = \sum_j \norm{d_j}^2 + \abs{a_0}^2. \]

%%
% One initilizes the algorithm with \(a_{J+1}=f\). The set of coefficients
% \(d_{j},a_j\) are computed from \(a_{j+1}\) as
% \[ \forall k=0,\ldots,2^j-1, \quad
%       a_j[k] = \frac{a_{j+1}[2k] + a_{j+1}[2k+1]}{\sqrt{2}}
%   \qandq 
%       d_j[k] = \frac{a_{j+1}[2k] - a_{j+1}[2k+1]}{\sqrt{2}}. \]


%%
% Shortcut to compute \(a_j, d_j\) from \(a_{j+1}\).


haar = @(a)[   a(1:2:length(a)) + a(2:2:length(a));
                    a(1:2:length(a)) - a(2:2:length(a)) ]/sqrt(2);   

%%
% Display the result of the one step of the transform.

clf;
subplot(2,1,1);
plot(f); axis('tight'); title('Signal');
subplot(2,1,2);
plot(haar(f)); axis('tight'); title('Transformed');

%%
% The output of the forward transform is stored in the variable |fw|.
% At a given iteration indexed by a scale |j|, it will store in |fw(1:2^(j+1))|
% the variable \(a_{j+1}\), and in the remaining entries the variable
% \(d_{j+1},d_{j+2},\ldots,d_J\).

%%
% Initializes the algorithm at scale \(j=J\).

j = J;

%%
% Initialize |fw| to the full signal. 

fw = f;

%%
% At iteration indexed by \(j\), 
% select the sub-part of the signal containing \(a_{j+1}\), 
% and apply it the Haar operator.

fw(1:2^(j+1)) = haar(fw(1:2^(j+1)));


%%
% Display the signal and its coarse coefficients.

s = 400; t = 40; 
clf;
subplot(2,1,1);
plot(f,'.-'); axis([s-t s+t 0 1]); title('Signal (zoom)');
subplot(2,1,2);
plot(fw(1:2^j),'.-'); axis([(s-t)/2 (s+t)/2 min(fw(1:2^j)) max(fw(1:2^j))]); title('Averages (zoom)');

%EXO
%% Implement a full wavelet Haar transform that extract iteratively wavelet
%% coefficients, by repeating these steps. Take care of choosing the
%% correct number of steps.
J = log2(N)-1;
Jmin = 0;
fw = f;
clf;
subplot(4,1,1);
plot(f); axis('tight'); title('Signal');
for j=J:-1:Jmin
    fw(1:2^(j+1)) = haar(fw(1:2^(j+1)));
    %    
    j1 = J-j;
    if j1<3
        d = fw(2^j+1:2^(j+1));
        subplot(4,1,j1+2);
        plot(1:2^(j1+1):N,d);  axis('tight');
        title( strcat(['Details, j=' num2str(j)]) );
    end   
end
%EXO

%% 
% Check that the transform is
% orthogonal, which means that the energy of the coefficient is the same
% as the energy of the signal.

fprintf('Should be 0: %.3f.\n', (norm(f)-norm(fw))/norm(f));

%%
% We display the whole set of coefficients |fw|, with red vertical
% separator between the scales.
% Can you recognize where are the low frequencies and the high frequencies? 

clf; 
plot_wavelet(fw);
axis([1 N -2 2]);

%% Backward Haar Transform
% The backward transform computes a signal \(f_1 = \Hh^*(h)\) from a set of
% coefficients \(  h = (d_j)_{j=0}^J \cup \{a_0\} \)

%%
% If \(h = \Hh(f)\) are the coefficients of a signal \(f\), one recovers
% \( f_1 = f \).

%%
% The inverse transform starts from \(j=0\), and computes \(a_{j+1}\)
% from the knowledge of \(a_j,d_j\) as
% \[
%   \forall k = 0,\ldots,2^j-1, \quad
%       a_{j+1}[2k] = \frac{ a_j[k] + d_j[k] }{\sqrt{2}}
%       a_{j+1}[2k+1] = \frac{ a_j[k] - d_j[k] }{\sqrt{2}}
% \]

%%
% One step of inverse transform

ihaar = @(a,d)assign( zeros(2*length(a),1), ...
        [1:2:2*length(a), 2:2:2*length(a)], [a+d; a-d]/sqrt(2) );


%%
% Initialize the signal to recover \(f_1\) as the transformed coefficient, and
% select the smallest possible scale \(j=0\).

f1 = fw;
j = 0;

%% 
% Apply one step of inverse transform.

f1(1:2^(j+1)) = ihaar(f1(1:2^j), f1(2^j+1:2^(j+1)));

%EXO
%% Write the inverse wavelet transform that computes |f1| from the
%% coefficients |fw|.
f1 = fw;
clf;
for j=Jmin:J
    f1(1:2^(j+1)) = ihaar(f1(1:2^j), f1(2^j+1:2^(j+1)));
    j1 = J-j;
    if j1<4
        subplot(4,1,j1+1);
        plot(1:2^j1:N,f1(1:2^(j+1)),'.-'); axis('tight');
        title( strcat(['Partial reconstruction, j=' num2str(j)]) );
    end
end
%EXO

%% 
% Check that we have correctly recovered the signal.

fprintf('Should be 0: %.3f.\n', (norm(f-f1))/norm(f));

%% Haar Linear Approximation
% Linear approximation is obtained by setting to zeros the
% Haar coefficient coefficients above a certain scale \(j\).

%%
% Cut-off scale.

j = J-1;

%%
% Set to zero fine scale coefficients.

fw1 = fw; fw1(2^j+1:end) = 0;

%%
% Computing the signal \(f_1\) corresponding to these coefficients |fw1|
% computes the orthogonal projection on the space of signal that are
% constant on disjoint intervals of length \(N/2^j\).

%EXO
%% Display the reconstructed signal obtained from |fw1|, for a decreasing cut-off scale \(j\).
jlist = J-(1:3);
fw = perform_haar_transf(f,1,+1);
clf;
for i=1:length(jlist)
    j = jlist(i);
    fw1 = fw; fw1(2^j+1:end) = 0;    
    f1 = perform_haar_transf(fw1,1,-1);
    % display
    subplot(length(jlist),1,i);
    hh = plot(f1); axis('tight');
    if using_matlab()
        set_linewidth(hh,2);
    end
    title( strcat(['j=' num2str(j) ', SNR=' num2str(snr(f,f1),3) 'dB']) );
end
%EXO


%% Haar Non-linear Approximation
% Non-linear approximation is obtained by thresholding low amplitude
% wavelet coefficients.
% It is defined as 
% \[ f_T = \Hh^* \circ S_T \circ \Hh(f) \]
% where \(S_T\) is the hard tresholding operator 
% \[ S_T(x)_i = \choice{
%       x_i \qifq \abs{x_i}>T, \\
%       0 \quad \text{otherwise}.
%   }. \]

S = @(x,T) x .* (abs(x)>T);

%%
% Set the threshold value.

T = .5;

%%
% Threshold the coefficients.

fwT = S(fw,T); 

%%
% Display the coefficients before and after thresholding.

clf;
subplot(2,1,1);
plot_wavelet(fw); axis('tight'); title('Original coefficients');
subplot(2,1,2);
plot_wavelet(fwT); axis('tight'); title('Thresholded coefficients');

%EXO
%% Find the threshold \(T\) so that the number of remaining coefficients in
%% |fwT| is a fixed number \(m\). Use this threshold to compute |fwT| and then display
%% the corresponding approximation \(f_1\) of \(f\). Try for an increasing number \(m\) of coeffiients.
% compute the threshold T
m_list = round([.05 .1 .2]*N); % number of kept coefficients
fw = perform_haar_transf(f,1,+1);
clf;
for i=1:length(m_list)
    m = m_list(i);
    % select threshold
    v = sort(abs(fw(:)));
    if v(1)<v(N)
        v = reverse(v);
    end
    T = v(m);
    fwT = fw .* (abs(fw)>=T);
    % inverse
    f1 = perform_haar_transf(fwT,1,-1);
    % display
    subplot(length(m_list),1,i);
    hh = plot(f1); axis('tight');
    if using_matlab()
        set_linewidth(hh,2);
    end
    title( strcat(['m=' num2str(m) ', SNR=' num2str(snr(f,f1),3) 'dB']) );
end
%EXO

%% The Shape of a Wavelet
% A wavelet coefficient corresponds to an inner product of \(f\) with a wavelet Haar atom
% \(\psi_{j,k}\)
% \[ d_j[k] = \dotp{f}{\psi_{j,k}} \]

%%
% The wavelet \(\psi_{j,k}\) can be computed by applying the inverse wavelet
% transform to |fw| where |fw[m]=1| and |fw[s]=0| for \(s \neq m\)
% where \(m = 2^j+k\).

%EXO
%% Compute wavelets at several positions and scales.
J = log2(N)-1;
selj = ( J-2:J )-3;
pos = [0 .5];
f = [];
clf;
k = 0;
for j=selj
    k = k+1;
    for q=1:length(pos)
        fw = zeros(N,1);
        p = 1 + (1+pos(q))*2^j;
        fw(p) = 1;
        f(:,q) = perform_haar_transf(fw,1,-1);
        f(:,q) = circshift(f(:,q),N/4);
    end
    f(1:N/2-1,2) = nan(); f(N/2+1:N,1) = nan();
    subplot(3,1,k);
    hh = plot(f); axis('tight');
    axis([1 N min(f(:))*1.05 max(f(:))*1.05]);
    if using_matlab()
        set_linewidth(hh,2);
    end
end
%EXO

%% Bibliography
% <html><a name="biblio"></a></html>

%%
% * [Haar1910] Haar A. <http://dx.doi.org/10.1007/BF01456927 _Zur Theorie der orthogonalen Funktionensysteme_>, Mathematische Annalen, 69, pp 331-371, 1910.

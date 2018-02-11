Jupyter Notebook
[M]ml_1_pca_nn
Last Checkpoint: 4 minutes ago
Autosave Failed!
Julia 0.5.0
Logout
Julia 0.5.0 No kernelerrorNot Trusted
File
Edit
View
Insert
Cell
Kernel
Help
Run
PCA, Nearest-Neighbors Classification and Clustering
Important: Please read the installation page for details about how to install the toolboxes.                                                       
This tour details Principal Component Analysis (dimentionality reduction), supervised classification using nearest neighbors and unsupervised clustering using kk-means.
We recommend that after doing this Numerical Tours, you apply it to your own data, for instance using a dataset from LibSVM.
Disclaimer: these machine learning tours are intended to be overly-simplistic implementations and applications of baseline machine learning methods. For more advanced uses and implementations, we recommend to use a state-of-the-art library, the most well known being Scikit-Learn.
In [497]:

1
using PyPlot
2
using NtToolBox
Dataset Loading
We use here the famous IRIS dataset of Fisher. The data set consists of 50 samples from each of three species of Iris (Iris setosa, Iris virginica and Iris versicolor). Four features were measured from each sample: the length and the width of the sepals and petals, in centimetres.
Helpers.
In [498]:

PlotBoxAspectRatio
1
# SetAR = @(ar)set(gca, 'PlotBoxAspectRatio', [1 ar 1], 'FontSize', 10);
Load the dataset.
In [499]:

1
A = readdlm("iris_dataset.csv", ',');
Randomly permute it.
In [500]:

1
A = A[randperm(size(A,1)),:];
In [501]:

','
1
A = readdlm("A", ',')
Out[501]:
150×5 Array{Float64,2}:
 7.7  2.8  6.7  2.0  2.0
 5.1  2.5  3.0  1.1  1.0
 5.4  3.4  1.5  0.4  0.0
 5.1  3.4  1.5  0.2  0.0
 5.1  3.7  1.5  0.4  0.0
 5.5  4.2  1.4  0.2  0.0
 6.1  3.0  4.6  1.4  1.0
 5.5  2.6  4.4  1.2  1.0
 6.7  3.0  5.2  2.3  2.0
 7.7  2.6  6.9  2.3  2.0
 6.4  2.7  5.3  1.9  2.0
 6.2  2.8  4.8  1.8  2.0
 4.9  3.1  1.5  0.1  0.0
 ⋮                      
 5.1  3.3  1.7  0.5  0.0
 5.6  2.7  4.2  1.3  1.0
 4.4  3.0  1.3  0.2  0.0
 4.8  3.0  1.4  0.1  0.0
 4.4  2.9  1.4  0.2  0.0
 6.7  3.1  4.4  1.4  1.0
 5.1  3.8  1.5  0.3  0.0
 6.3  3.3  4.7  1.6  1.0
 5.6  2.8  4.9  2.0  2.0
 4.9  3.1  1.5  0.1  0.0
 4.8  3.4  1.6  0.2  0.0
 7.7  3.8  6.7  2.2  2.0
Separate the features (xi)ni=1(xi)i=1nfrom the class information. The feature are stored as the row of a matrix X∈ℝn×pX∈Rn×p Be sure to start the class at index 1.
In [502]:

)
1
X = A[:,1:end-1]
2
y = Int.(A[:,end])
3
y = y-minimum(y)+1;
nn is the number of samples, pp is the dimensionality of the features, kk is the number of classes.
In [503]:

1
n,p = size(X)
2
k = maximum(y);
Dimensionality Reduction and PCA
In order to display in 2-D or 3-D the data, dimensionality reduction is needed. The simplest method is the Principal Component Analysis (PCA), which perform an orthogonal linear projection on the principal axsis (eigenvector) of the covariance matrix.
Compute empirical mean
m=1n∑i=1nxi∈ℝp
m=1n∑i=1nxi∈Rp
and covariance
C=1n∑i=1n(xi−m)(xi−m)⊤∈ℝp×p.
C=1n∑i=1n(xi−m)(xi−m)⊤∈Rp×p.
Denoting X̃ =X−1pm⊤X~=X−1pm⊤, one has C=X̃ ⊤X̃ C=X~⊤X~.
In [504]:

1
Xm = X -> X-repeat(mean(X,1), outer=(size(X,1), 1))
2
Cov = X -> Xm(X)'*Xm(X);
Display the covariance matrix.
In [505]:

cmap = get_cmap("jet")
1
imshow(Cov(X), extent=[0, 1, 0, 1], cmap = get_cmap("jet"));

Compute PCA ortho-basis using the SVD decomposition
X̃ =Udiag(d)V
X~=Udiag(d)V
where U∈ℝn×pU∈Rn×p and V∈ℝp×pV∈Rp×p have orthonormal columns. VV are the principal directions of variance and are order by decreasing variances.
In [506]:

1
U,D,V = svd(Xm(X),thin=true);
Compute the feature in the PCA basis, zi=V⊤(xi−m)zi=V⊤(xi−m), stored in matrix format as Z=X̃ VZ=X~V.
In [507]:

1
Z = Xm(X) * V;
Plot the singular values of the covariances, which corresponds to the standard deviation of the data along the principal directions.
In [508]:

1
#figure(figsize = (5,3))
2
plot(D, ".-", linewidth= 2, markersize= 20);

The first dimensions of the zizi are the optimal way to linearly embed the data in a low dimensional space. This can be used for display in 2-D using the first two dimension.
In [509]:

1
figure(figsize=(7,3))
2
col = [ [1 0 0]; [0 1 0]; [0 0 1]; [0 0 0]; [0 1 1]; [1 0 1]; [1 1 0]; 
3
    [1 .5 .5]; [.5 1 .5]; [.5 .5 1]  ]';
4
ms = 25;
5
lgd = [];
6
for i=1:min(k,size(col,2))
7
    I = find(y.==i)
8
    plot(Z[I,1], Z[I,2], ".", markersize= ms, c=col[:,Int(i)], markersize=ms);
9
    append!(lgd,i);
10
end
11
#axis tight; axis equal; box on;
12
#legend(lgd, "Location", "EastOutside");
13
#SetAR(1);

Similar display in 3-D.
In [510]:

1
#figure(figsize=(10,5))
2
for i=1:k
3
    I = find(y.==i)
4
    plot3D(Z[I,1], Z[I,2], Z[I,3], ".", c=col[:,Int(i)], markersize=ms, label=Int(i))
5
end
6
gca()[:view_init](50, 250)
7
axis("equal"); axis("tight"); gca()[:grid](false);
8
legend(loc="best");

Supervised Learning: Nearest Neighbor Classification
Probably the simplest method for supervised classification is Nearest Neighbor (RR-NN), where RR is a parameter indexing the number of neighbor. Increasing RR is important to cope with noise and obtain smoother decision boundary, and hence better generalization performance.
The class predicted for a point xx is the one which is the most represented among the RR points (xi)i(xi)i which are the closed to xx.
Split into training and testing.
In [511]:

1
n0 = Int(round(.5*n)); n1 = n-n0
2
X0 = X[1:n0,:];     y0 = y[1:n0]
3
X1 = X[n0+1:end,:]; y1 = y[n0+1:end];
Macro to compute pairwise squared Euclidean distance matrix.
In [512]:

-2*(X*Z')
1
distmat = (X,Z) -> broadcast(+,sum(X'.*X',1)',sum(Z'.*Z',1))-2*(X*Z');
Compute Euclidean distance between some xx and all other x1,jx1,j in the training set.
In [513]:

D = distmat(X0,x);
1
i = 1; x = X1[i,:]';  # could be any point
2
D = distmat(X0,x);
Sort the distance and generate the list of sorted classes yσ=(yσ(i))iyσ=(yσ(i))i. This generate an indexing σσ (a permutation of {1,…,n}{1,…,n}) such that
∥x−xσ(1)∥≤∥x−xσ(2)∥≤…≤∥x−xσ(n)∥.
∥x−xσ(1)∥≤∥x−xσ(2)∥≤…≤∥x−xσ(n)∥.
In [514]:

1
I = sortperm(D[:])
2
ys = Int.(y[I]);
For a given RR, one can compute the histogram of class apparition
hℓ≡1R{i,σ(i)∈{1,…,R}}=♯σ−1({1,…,R}).
hℓ≡1R{i,σ(i)∈{1,…,R}}=♯σ−1({1,…,R}).
The decision class for xx is then the maximum of the histogram
c(x)≡argmaxℓhℓ
c(x)≡argmaxℓhℓ
In [515]:

1
function custom_hist(h)
2
    return [sum(h .== i) for i in 1:k]
3
end
WARNING: Method definition custom_hist(Any) in module Main at In[491]:2 overwritten at In[515]:2.
Out[515]:
custom_hist (generic function with 1 method)
In [516]:

1
R = 5
2
h = custom_hist(ys[1:R,:])
3
h = h / R
4
c = indmax(h);
5
print("c(x)=",c," [true class=",Int(y1[i]), "]\n");
c(x)=3 [true class=3]
Display the histigram (hℓ)ℓ(hℓ)ℓ of reparttion of class indexes as RR grows.
In [517]:

])
1
figure(figsize=(8,6))
2
Rlist = round([.05 .1 .5 1]*n0); #  [5 50 100]
3
clf()
4
for i=1:length(Rlist)
5
    R = Int(Rlist[i]);
6
    h = custom_hist(ys[1:R,:])
7
    h = h / R;
8
    subplot(length(Rlist),1,i);
9
    bar(1:k,h); 
10
    axis([0.5 k+.5 0 1]');
11
end

Exercise 1
Perform the NN classification for all the points in the test set, and for varying RR. Show how the classification score SS (number of correctly classified points) evolves with RR lot(1:Rmax, S, '.-', 'MarkerSize', ms);
In [520]:

1
1
include("NtSolutions/ml_1_pca_nn/exo1.jl");

In [521]:

#
1
# Insert your code here.
Exercise 2
Display, as a function of the position in 2-D PCA space, the class output by the RR-NN method when applied in 2-D. ounding boxes est for different R
In [534]:

1
?meshgrid
search: meshgrid

Out[534]:
No documentation found.
NtToolBox.meshgrid is a Function.
# 3 methods for generic function "meshgrid":
meshgrid(v::AbstractArray{T<:Any,1}) at /Users/quentin/.julia/v0.5/NtToolBox/src/ndgrid.jl:33
meshgrid{T}(vx::AbstractArray{T,1}, vy::AbstractArray{T,1}) at /Users/quentin/.julia/v0.5/NtToolBox/src/ndgrid.jl:36
meshgrid{T}(vx::AbstractArray{T,1}, vy::AbstractArray{T,1}, vz::AbstractArray{T,1}) at /Users/quentin/.julia/v0.5/NtToolBox/src/ndgrid.jl:44
In [538]:

; Z = z1;
1
X=Z[:,1:2]; Z = z1;
In [540]:

1
broadcast(+,sum(X'.*X',1)',sum(Z'.*Z',1))-2*(X*Z')
Out[540]:
150×40000 Array{Float64,2}:
 11.2289   11.2079   11.1897   …  69.6207  70.1698  70.7218  71.2767
 42.7881   42.4308   42.0764      16.9435  17.1564  17.3722  17.5909
 49.9006   49.4288   48.9599      19.3608  19.4591  19.5604  19.6645
 53.4322   52.9466   52.464       17.4349  17.5195  17.607   17.6974
 51.423    50.941    50.4619      19.2257  19.3139  19.405   19.499 
 48.1169   47.6308   47.1477   …  25.1313  25.2154  25.3024  25.3922
 23.0937   22.8736   22.6564      35.033   35.3831  35.7361  36.092 
 31.0291   30.7765   30.5267      26.9336  27.2511  27.5715  27.8948
 16.4374   16.2977   16.1609      47.4898  47.9202  48.3536  48.7899
 12.5417   12.5431   12.5475      72.261   72.8326  73.4071  73.9845
 20.0563   19.9058   19.7582   …  42.865   43.2846  43.7071  44.1326
 22.228    22.0359   21.8467      37.4551  37.8331  38.2141  38.5981
 57.0627   56.5707   56.0817      14.6142  14.6924  14.7736  14.8576
  ⋮                            ⋱                                    
 50.7933   50.3296   49.8689      17.1082  17.2147  17.3241  17.4364
 30.2998   30.0389   29.7809      27.1215  27.4308  27.743   28.0581
 64.1964   63.681    63.1685   …  11.4733  11.528   11.5857  11.6462
 59.5042   59.0036   58.5059      13.5203  13.5898  13.6623  13.7377
 63.7152   63.2069   62.7016      11.0212  11.0831  11.1479  11.2156
 19.2856   19.0684   18.8541      39.4048  39.7577  40.1135  40.4723
 51.4508   50.9655   50.4831      19.8151  19.8999  19.9876  20.0783
 19.6128   19.4083   19.2068   …  39.4402  39.8059  40.1745  40.546 
 26.1046   25.908    25.7143      33.9689  34.3425  34.7189  35.0983
 57.0627   56.5707   56.0817      14.6142  14.6924  14.7736  14.8576
 55.2971   54.8098   54.3254      15.7022  15.7851  15.8708  15.9595
  6.97446   6.95258   6.93361     76.0586  76.6069  77.158   77.7121
In [*]:

R
1
# bounding boxes
2
B = maximum(abs(Z[:,1:2]))
3
q = 200
4
r = linspace(-B,B,q)
5
V,U = meshgrid(r,r)
6
z1 = [U[:] V[:]]
7
#test for different R
8
Rlist = [1 5 10 40]
9
for ir=1:length(Rlist)
10
    R=Rlist[ir]
11
    print(R)
12
    #
13
    D = distmat(Z[:,1:2],z1);
14
    Ds = sort(D,1)
15
    I = mapslices(sortperm, D, 1)
16
    ys = y[I]
17
    #
18
    if R==1
19
        C = ys[1,:]
20
    else
21
        h = mapslices(custom_hist, ys[1:R,:],1)
22
        C = mapslices(indmax, h,1)
23
    end
24
    C = reshape(C, outer=(q,q))
25
    # maps class to color
26
    Cr = zeros(q,q,3)
27
    for i=1:k
28
        for a=1:3
29
            Cr[:,:,a] = Cr[:,:,a] + (C.==i)*col(a,i);
30
        end
31
    end
32
    # display
33
    subplot(2,2,ir)
34
    imagesc(r,r,permute(Cr,[2 1 3]), cmap = get_cmap("jet"));
35
    for i=1:k
36
        I = find(y==i);
37
        plot(Z[I,1], Z[I,2], "o", c=col(:,i)*.9, cmap = get_cmap("jet"))#, 'MarkerSize', 5, 'MarkerEdgeColor', col(:,i)*.5);
38
    #    plot(Z(I,1), Z(I,2), '.', 'Color', col(:,i), 'MarkerSize', ms);
39
    end
40
    axis("tight"); axis("equal"); axis("off")
41
    title(["R=" R]);
42
end
43
​
In [22]:

1
exo2()

In [23]:

1
%% Insert your code here.
Unsupervised Learning: kk-means
In an un-supervised setting, the class information yy is not available. The basic problem is then to recover class information from the knowledge of xx only. This corresponds to the clustering problem.
Select a subset of classes
In [24]:

1
if k>=4
2
ksvg = k; Xsvg = X; ysvg = y;
3
k = 3;
4
I = find(y<=k);
5
X = X(I,:); y = y(I);
6
n = length(I);
7
end
PCA
In [25]:

1
[U,D,V] = svd(Xm(X),'econ');
2
Z = Xm(X) * V;
The most basic algorithm is the , which tries to recover the class index y¯i=ℓy¯i=ℓ from the distance ∥xi−cℓ∥∥xi−cℓ∥ between the feature point xixi and the class centroid cℓcℓ (which are the unknown of the problem).
It does so by minimizing the following non-convex energy
min(cℓ)ℓ∑iminℓ∥xi−cℓ∥2
min(cℓ)ℓ∑iminℓ∥xi−cℓ∥2
We first initialize the class centroids (cℓ)ℓ(cℓ)ℓ at random among the points. They are stored in as the row of a matrix C∈ℝk×pC∈Rk×p.
In [26]:

1
I = randperm(n); I = I(1:k);
2
C = X(I,:);
The kk-means algorithm iterate between first determining the class of each point using the distance to the centroids
∀i∈{1,…,n},y¯i←argminℓ∥xi−cℓ∥.
∀i∈{1,…,n},y¯i←argminℓ∥xi−cℓ∥.
In [27]:

1
D = distmat(X,C);
2
[~,yb] = min(D, [], 2);
Display the centroids and the classes using colors. This correspodns to a Voronoi diagram segmentation in the high dimensional space, but here the display is done in 2D.
In [28]:

1
clf;
2
hold on;
3
for i=1:k
4
    I = find(yb==i);
5
    plot(Z(I,1), Z(I,2), '.', 'Color', col(:,i), 'MarkerSize', 25);
6
end
7
CV = (C-repmat(mean(X,1), [k 1]))*V;
8
for i=1:k
9
    plot(CV(i,1), CV(i,2), 'o', 'MarkerFaceColor', col(:,i), 'MarkerSize', 12, 'MarkerEdgeColor', 'k');
10
end
11
axis tight; axis equal; axis off;
12
SetAR(1);

The second step of the kk-means algorithm is to update the centroids position to be the mean of the points inside each class
∀ℓ∈{1,…,k},cℓ←∑i:yi=ℓxi♯{i:yi=ℓ}.
∀ℓ∈{1,…,k},cℓ←∑i:yi=ℓxi♯{i:yi=ℓ}.
In [29]:

1
for l=1:k
2
    C(l,:) = mean( X(yb==l,:), 1 );
3
end
Exercise 3
Peform several step of the kk-means algorithm. nit
In [30]:

1
exo3()

In [31]:

1
%% Insert your code here.
Display the histogram of (true, i.e. according to yy) class inside each estimated class (i.e. according to y¯y¯).
In [32]:

1
clf
2
for l=1:k
3
    I = find(yb==l);
4
    h = hist(y(I),1:k); h = h/sum(h);
5
    subplot(k,1,l);
6
    bar(1:k,h); 
7
    axis([0.5 k+.5 0 1]);
8
    set(gca, 'FontSize', 10);
9
end

Exercise 4
Implement better initialization strategies such as farthest point sampling or .
In [33]:

1
exo4()
In [34]:

1
%% Insert your code here.

1
<script>
2
  $(document).ready(function(){
3
      $('div.prompt').hide();
4
  });
5
</script>

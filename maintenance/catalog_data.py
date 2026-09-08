"""Editorial metadata shared by notebooks and the website catalogue."""

REFERENCES = {
    "mallat": (
        "Stéphane Mallat",
        "A Wavelet Tour of Signal Processing: The Sparse Way",
        "2009, 3rd ed., Academic Press",
        "https://www.di.ens.fr/~mallat/biblio.html",
        "Multiresolution analysis, sparse approximation, and wavelet algorithms.",
    ),
    "daubechies": (
        "Ingrid Daubechies",
        "Ten Lectures on Wavelets",
        "1992, SIAM",
        "https://doi.org/10.1137/1.9781611970104",
        "Compactly supported orthogonal wavelets and their regularity.",
    ),
    "peyre": (
        "Gabriel Peyré",
        "Advanced Signal, Image and Surface Processing",
        "2010, course notes",
        "https://www.numerical-tours.com/book/AdvancedSignalProcessing.pdf",
        "A mathematical companion to the Numerical Tours.",
    ),
    "strang": (
        "Gilbert Strang",
        "Linear Algebra and Learning from Data",
        "2019, Wellesley-Cambridge Press",
        "https://math.mit.edu/~gs/learningfromdata/",
        "Matrix factorizations, least squares, and low-rank representations.",
    ),
    "scipy": (
        "Pauli Virtanen et al.",
        "SciPy 1.0: Fundamental Algorithms for Scientific Computing in Python",
        "2020, Nature Methods 17, 261–272",
        "https://doi.org/10.1038/s41592-019-0686-2",
        "The numerical routines used for transforms, interpolation, and optimization.",
    ),
    "skimage": (
        "Stéfan van der Walt et al.",
        "scikit-image: image processing in Python",
        "2014, PeerJ 2:e453",
        "https://doi.org/10.7717/peerj.453",
        "Reproducible image processing and image measurements in Python.",
    ),
    "boyd": (
        "Stephen Boyd and Lieven Vandenberghe",
        "Convex Optimization",
        "2004, Cambridge University Press",
        "https://web.stanford.edu/~boyd/cvxbook/",
        "Convexity, duality, optimality conditions, and interior-point methods.",
    ),
    "prox": (
        "Neal Parikh and Stephen Boyd",
        "Proximal Algorithms",
        "2014, Foundations and Trends in Optimization 1(3), 127–239",
        "https://web.stanford.edu/~boyd/papers/prox_algs.html",
        "Proximity operators and splitting methods for nonsmooth objectives.",
    ),
    "fista": (
        "Amir Beck and Marc Teboulle",
        "A Fast Iterative Shrinkage-Thresholding Algorithm for Linear Inverse Problems",
        "2009, SIAM Journal on Imaging Sciences 2(1), 183–202",
        "https://doi.org/10.1137/080716542",
        "Accelerated proximal gradient descent and its convergence rate.",
    ),
    "chambollepock": (
        "Antonin Chambolle and Thomas Pock",
        "A First-Order Primal-Dual Algorithm for Convex Problems with Applications to Imaging",
        "2011, Journal of Mathematical Imaging and Vision 40, 120–145",
        "https://doi.org/10.1007/s10851-010-0251-1",
        "Primal–dual splitting for total variation and other composite penalties.",
    ),
    "condat": (
        "Laurent Condat",
        "A Primal–Dual Splitting Method for Convex Optimization Involving Lipschitzian, Proximable and Linear Composite Terms",
        "2013, Journal of Optimization Theory and Applications 158, 460–479",
        "https://doi.org/10.1007/s10957-012-0245-9",
        "Step-size conditions and a general primal–dual framework.",
    ),
    "rof": (
        "Leonid Rudin, Stanley Osher, and Emad Fatemi",
        "Nonlinear Total Variation Based Noise Removal Algorithms",
        "1992, Physica D 60, 259–268",
        "https://doi.org/10.1016/0167-2789(92)90242-F",
        "The foundational total-variation image restoration model.",
    ),
    "donoho": (
        "David L. Donoho",
        "De-noising by Soft-Thresholding",
        "1995, IEEE Transactions on Information Theory 41(3), 613–627",
        "https://doi.org/10.1109/18.382009",
        "Why shrinkage of wavelet coefficients suppresses noise.",
    ),
    "stein": (
        "Charles M. Stein",
        "Estimation of the Mean of a Multivariate Normal Distribution",
        "1981, Annals of Statistics 9(6), 1135–1151",
        "https://doi.org/10.1214/aos/1176345632",
        "Unbiased risk estimation for Gaussian observations.",
    ),
    "cai": (
        "T. Tony Cai",
        "Adaptive Wavelet Estimation: A Block Thresholding and Oracle Inequality Approach",
        "1999, Annals of Statistics 27(3), 898–924",
        "https://doi.org/10.1214/aos/1018031262",
        "Statistical motivation for processing neighboring coefficients together.",
    ),
    "nlm": (
        "Antoni Buades, Bartomeu Coll, and Jean-Michel Morel",
        "A Non-Local Algorithm for Image Denoising",
        "2005, CVPR, vol. 2, 60–65",
        "https://doi.org/10.1109/CVPR.2005.38",
        "Patch similarity as the basis for nonlocal averaging.",
    ),
    "nlmreview": (
        "Antoni Buades, Bartomeu Coll, and Jean-Michel Morel",
        "A Review of Image Denoising Algorithms, with a New One",
        "2005, Multiscale Modeling & Simulation 4(2), 490–530",
        "https://doi.org/10.1137/040616024",
        "Comparison of local, transform, and nonlocal denoising models.",
    ),
    "cover": (
        "Thomas M. Cover and Joy A. Thomas",
        "Elements of Information Theory",
        "2006, 2nd ed., Wiley",
        "https://doi.org/10.1002/047174882X",
        "Entropy, source coding, and fundamental compression limits.",
    ),
    "shannon": (
        "Claude E. Shannon",
        "A Mathematical Theory of Communication",
        "1948, Bell System Technical Journal 27, 379–423 and 623–656",
        "https://doi.org/10.1002/j.1538-7305.1948.tb01338.x",
        "Entropy and the source coding theorem.",
    ),
    "huffman": (
        "David A. Huffman",
        "A Method for the Construction of Minimum-Redundancy Codes",
        "1952, Proceedings of the IRE 40(9), 1098–1101",
        "https://doi.org/10.1109/JRPROC.1952.273898",
        "The greedy construction of optimal binary prefix codes.",
    ),
    "esl": (
        "Trevor Hastie, Robert Tibshirani, and Jerome Friedman",
        "The Elements of Statistical Learning",
        "2009, 2nd ed., Springer",
        "https://hastie.su.domains/ElemStatLearn/",
        "Regression, classification, regularization, and model assessment.",
    ),
    "lasso": (
        "Robert Tibshirani",
        "Regression Shrinkage and Selection via the Lasso",
        "1996, Journal of the Royal Statistical Society B 58(1), 267–288",
        "https://doi.org/10.1111/j.2517-6161.1996.tb02080.x",
        "Sparse regression through an absolute-value penalty.",
    ),
    "sklearn": (
        "Fabian Pedregosa et al.",
        "Scikit-learn: Machine Learning in Python",
        "2011, Journal of Machine Learning Research 12, 2825–2830",
        "https://www.jmlr.org/papers/v12/pedregosa11a.html",
        "Practical estimators and reproducible model evaluation.",
    ),
    "deep": (
        "Ian Goodfellow, Yoshua Bengio, and Aaron Courville",
        "Deep Learning",
        "2016, MIT Press",
        "https://www.deeplearningbook.org/",
        "Neural networks, automatic differentiation, and optimization.",
    ),
    "adam": (
        "Diederik P. Kingma and Jimmy Ba",
        "Adam: A Method for Stochastic Optimization",
        "2015, ICLR",
        "https://arxiv.org/abs/1412.6980",
        "Adaptive stochastic gradient updates for neural network training.",
    ),
    "backprop": (
        "David E. Rumelhart, Geoffrey E. Hinton, and Ronald J. Williams",
        "Learning Representations by Back-Propagating Errors",
        "1986, Nature 323, 533–536",
        "https://doi.org/10.1038/323533a0",
        "The chain-rule computation underlying multilayer network training.",
    ),
    "sgd": (
        "Léon Bottou, Frank E. Curtis, and Jorge Nocedal",
        "Optimization Methods for Large-Scale Machine Learning",
        "2018, SIAM Review 60(2), 223–311",
        "https://arxiv.org/abs/1606.04838",
        "Stochastic gradients, variance reduction, and practical convergence.",
    ),
    "ot": (
        "Gabriel Peyré and Marco Cuturi",
        "Computational Optimal Transport",
        "2019, Foundations and Trends in Machine Learning 11(5–6), 355–607",
        "https://arxiv.org/abs/1803.00567",
        "Transport plans, duality, entropic regularization, and applications.",
    ),
    "sinkhorn": (
        "Marco Cuturi",
        "Sinkhorn Distances: Lightspeed Computation of Optimal Transport",
        "2013, NeurIPS",
        "https://arxiv.org/abs/1306.0895",
        "Entropic regularization and scalable matrix scaling.",
    ),
    "scaling": (
        "Richard Sinkhorn and Paul Knopp",
        "Concerning Nonnegative Matrices and Doubly Stochastic Matrices",
        "1967, Pacific Journal of Mathematics 21(2), 343–348",
        "https://doi.org/10.2140/pjm.1967.21.343",
        "The convergence foundations of alternating matrix normalization.",
    ),
    "unbalanced": (
        "Lénaïc Chizat, Gabriel Peyré, Bernhard Schmitzer, and François-Xavier Vialard",
        "Scaling Algorithms for Unbalanced Transport Problems",
        "2018, Mathematics of Computation 87, 2563–2609",
        "https://arxiv.org/abs/1607.05816",
        "Relaxed marginal constraints and generalized Sinkhorn iterations.",
    ),
    "semidiscrete": (
        "Quentin Mérigot",
        "A Multiscale Approach to Optimal Transport",
        "2011, Computer Graphics Forum 30(5), 1583–1592",
        "https://doi.org/10.1111/j.1467-8659.2011.02032.x",
        "Geometric algorithms for semi-discrete transport.",
    ),
    "benamou": (
        "Jean-David Benamou, Guillaume Carlier, Marco Cuturi, Luca Nenna, and Gabriel Peyré",
        "Iterative Bregman Projections for Regularized Transportation Problems",
        "2015, SIAM Journal on Scientific Computing 37(2), A1111–A1138",
        "https://arxiv.org/abs/1412.5154",
        "Projection-based algorithms for transport and barycenters.",
    ),
    "canny": (
        "John Canny",
        "A Computational Approach to Edge Detection",
        "1986, IEEE Transactions on Pattern Analysis and Machine Intelligence 8(6), 679–698",
        "https://doi.org/10.1109/TPAMI.1986.4767851",
        "Detection, localization, and nonmaximum suppression of edges.",
    ),
    "snakes": (
        "Michael Kass, Andrew Witkin, and Demetri Terzopoulos",
        "Snakes: Active Contour Models",
        "1988, International Journal of Computer Vision 1, 321–331",
        "https://doi.org/10.1007/BF00133570",
        "Parametric curves driven by regularization and image forces.",
    ),
    "levelset": (
        "Stanley Osher and James A. Sethian",
        "Fronts Propagating with Curvature-Dependent Speed: Algorithms Based on Hamilton–Jacobi Formulations",
        "1988, Journal of Computational Physics 79(1), 12–49",
        "https://doi.org/10.1016/0021-9991(88)90002-2",
        "Implicit front propagation using level-set functions.",
    ),
    "chanvese": (
        "Tony F. Chan and Luminita A. Vese",
        "Active Contours Without Edges",
        "2001, IEEE Transactions on Image Processing 10(2), 266–277",
        "https://doi.org/10.1109/83.902291",
        "Region-based segmentation when image gradients are weak.",
    ),
    "fastmarch": (
        "James A. Sethian",
        "A Fast Marching Level Set Method for Monotonically Advancing Fronts",
        "1996, PNAS 93(4), 1591–1595",
        "https://doi.org/10.1073/pnas.93.4.1591",
        "A causal numerical solver for the eikonal equation.",
    ),
    "dijkstra": (
        "Edsger W. Dijkstra",
        "A Note on Two Problems in Connexion with Graphs",
        "1959, Numerische Mathematik 1, 269–271",
        "https://doi.org/10.1007/BF01386390",
        "Shortest paths with nonnegative edge weights.",
    ),
    "isomap": (
        "Joshua B. Tenenbaum, Vin de Silva, and John C. Langford",
        "A Global Geometric Framework for Nonlinear Dimensionality Reduction",
        "2000, Science 290, 2319–2323",
        "https://doi.org/10.1126/science.290.5500.2319",
        "Geodesic distances and multidimensional scaling for manifold learning.",
    ),
    "geometry": (
        "Mario Botsch, Leif Kobbelt, Mark Pauly, Pierre Alliez, and Bruno Lévy",
        "Polygon Mesh Processing",
        "2010, AK Peters",
        "https://www.pmp-book.org/",
        "Discrete differential operators, smoothing, and parameterization.",
    ),
    "floater": (
        "Michael S. Floater",
        "Parametrization and Smooth Approximation of Surface Triangulations",
        "1997, Computer Aided Geometric Design 14(3), 231–250",
        "https://doi.org/10.1016/S0167-8396(96)00031-3",
        "Barycentric mappings of a mesh into a planar domain.",
    ),
    "desbrun": (
        "Mathieu Desbrun, Mark Meyer, Peter Schröder, and Alan H. Barr",
        "Implicit Fairing of Irregular Meshes Using Diffusion and Curvature Flow",
        "1999, SIGGRAPH, 317–324",
        "https://doi.org/10.1145/311535.311576",
        "Laplacian smoothing and stable geometric diffusion.",
    ),
    "stam": (
        "Jos Stam",
        "Stable Fluids",
        "1999, SIGGRAPH, 121–128",
        "https://doi.org/10.1145/311535.311548",
        "Semi-Lagrangian advection and incompressible velocity projection.",
    ),
    "bridson": (
        "Robert Bridson",
        "Fluid Simulation for Computer Graphics",
        "2015, 2nd ed., CRC Press",
        "https://doi.org/10.1201/b18320",
        "Discretization, advection, pressure solves, and boundary conditions.",
    ),
    "galerne": (
        "Bruno Galerne, Yann Gousseau, and Jean-Michel Morel",
        "Random Phase Textures: Theory and Synthesis",
        "2011, IEEE Transactions on Image Processing 20(1), 257–267",
        "https://doi.org/10.1109/TIP.2010.2052822",
        "Stationary Gaussian textures and Fourier phase randomization.",
    ),
    "gatys": (
        "Leon A. Gatys, Alexander S. Ecker, and Matthias Bethge",
        "Texture Synthesis Using Convolutional Neural Networks",
        "2015, NeurIPS",
        "https://arxiv.org/abs/1505.07376",
        "Texture statistics from Gram matrices of deep features.",
    ),
    "portilla": (
        "Javier Portilla and Eero P. Simoncelli",
        "A Parametric Texture Model Based on Joint Statistics of Complex Wavelet Coefficients",
        "2000, International Journal of Computer Vision 40, 49–70",
        "https://doi.org/10.1023/A:1026553619983",
        "A multiscale statistical model for texture synthesis.",
    ),
    "ddpm": (
        "Jonathan Ho, Ajay Jain, and Pieter Abbeel",
        "Denoising Diffusion Probabilistic Models",
        "2020, NeurIPS",
        "https://arxiv.org/abs/2006.11239",
        "Learning a generative reverse process by denoising.",
    ),
    "scoresde": (
        "Yang Song et al.",
        "Score-Based Generative Modeling through Stochastic Differential Equations",
        "2021, ICLR",
        "https://arxiv.org/abs/2011.13456",
        "Reverse-time SDEs and probability-flow ODEs.",
    ),
    "vincent": (
        "Pascal Vincent",
        "A Connection Between Score Matching and Denoising Autoencoders",
        "2011, Neural Computation 23(7), 1661–1674",
        "https://doi.org/10.1162/NECO_a_00142",
        "The denoising score-matching objective used in diffusion training.",
    ),
    "hyvarinen": (
        "Aapo Hyvärinen",
        "Estimation of Non-Normalized Statistical Models by Score Matching",
        "2005, JMLR 6, 695–709",
        "https://www.jmlr.org/papers/v6/hyvarinen05a.html",
        "Estimating log-density gradients without a normalization constant.",
    ),
    "conformal": (
        "Anastasios N. Angelopoulos and Stephen Bates",
        "A Gentle Introduction to Conformal Prediction and Distribution-Free Uncertainty Quantification",
        "2023, Foundations and Trends in Machine Learning 16(4), 494–591",
        "https://arxiv.org/abs/2107.07511",
        "Coverage guarantees, exchangeability, and practical calibration.",
    ),
    "vovk": (
        "Vladimir Vovk, Alex Gammerman, and Glenn Shafer",
        "Algorithmic Learning in a Random World",
        "2005, Springer",
        "https://doi.org/10.1007/b106715",
        "The foundational theory of conformal prediction.",
    ),
    "shafer": (
        "Glenn Shafer and Vladimir Vovk",
        "A Tutorial on Conformal Prediction",
        "2008, JMLR 9, 371–421",
        "https://www.jmlr.org/papers/v9/shafer08a.html",
        "Conformity scores, p-values, and prediction regions.",
    ),
    "lei": (
        "Jing Lei et al.",
        "Distribution-Free Predictive Inference for Regression",
        "2018, JASA 113(523), 1094–1111",
        "https://arxiv.org/abs/1604.04173",
        "Regression prediction bands and split conformal inference.",
    ),
    "cqr": (
        "Yaniv Romano, Evan Patterson, and Emmanuel Candès",
        "Conformalized Quantile Regression",
        "2019, NeurIPS",
        "https://arxiv.org/abs/1905.03222",
        "Adaptive prediction intervals for heteroscedastic regression.",
    ),
    "graphlasso": (
        "Jerome Friedman, Trevor Hastie, and Robert Tibshirani",
        "Sparse Inverse Covariance Estimation with the Graphical Lasso",
        "2008, Biostatistics 9(3), 432–441",
        "https://doi.org/10.1093/biostatistics/kxm045",
        "Coordinate descent for sparse Gaussian precision matrices.",
    ),
    "banerjee": (
        "Onureena Banerjee, Laurent El Ghaoui, and Alexandre d’Aspremont",
        "Model Selection Through Sparse Maximum Likelihood Estimation for Multivariate Gaussian or Binary Data",
        "2008, JMLR 9, 485–516",
        "https://www.jmlr.org/papers/v9/banerjee08a.html",
        "The convex sparse covariance estimation problem.",
    ),
    "meinshausen": (
        "Nicolai Meinshausen and Peter Bühlmann",
        "High-Dimensional Graphs and Variable Selection with the Lasso",
        "2006, Annals of Statistics 34(3), 1436–1462",
        "https://doi.org/10.1214/009053606000000281",
        "Neighborhood selection and conditional-independence graphs.",
    ),
    "completion": (
        "Emmanuel J. Candès and Benjamin Recht",
        "Exact Matrix Completion via Convex Optimization",
        "2009, Foundations of Computational Mathematics 9, 717–772",
        "https://arxiv.org/abs/0805.4471",
        "Low-rank recovery through nuclear-norm minimization.",
    ),
    "svt": (
        "Jian-Feng Cai, Emmanuel J. Candès, and Zuowei Shen",
        "A Singular Value Thresholding Algorithm for Matrix Completion",
        "2010, SIAM Journal on Optimization 20(4), 1956–1982",
        "https://arxiv.org/abs/0810.3286",
        "The proximal operator of the nuclear norm.",
    ),
    "cs": (
        "David L. Donoho",
        "Compressed Sensing",
        "2006, IEEE Transactions on Information Theory 52(4), 1289–1306",
        "https://doi.org/10.1109/TIT.2006.871582",
        "Sparse recovery from undersampled linear measurements.",
    ),
    "ksvd": (
        "Michal Aharon, Michael Elad, and Alfred Bruckstein",
        "K-SVD: An Algorithm for Designing Overcomplete Dictionaries for Sparse Representation",
        "2006, IEEE Transactions on Signal Processing 54(11), 4311–4322",
        "https://doi.org/10.1109/TSP.2006.881199",
        "Alternating sparse coding and dictionary updates.",
    ),
    "mairal": (
        "Julien Mairal, Francis Bach, Jean Ponce, and Guillermo Sapiro",
        "Online Learning for Matrix Factorization and Sparse Coding",
        "2010, JMLR 11, 19–60",
        "https://www.jmlr.org/papers/v11/mairal10a.html",
        "Scalable dictionary learning and matrix factorization.",
    ),
    "bss": (
        "Özgür Yılmaz and Scott Rickard",
        "Blind Separation of Speech Mixtures via Time-Frequency Masking",
        "2004, IEEE Transactions on Signal Processing 52(7), 1830–1847",
        "https://doi.org/10.1109/TSP.2004.828896",
        "Separating audio sources through sparse time-frequency structure.",
    ),
    "ica": (
        "Aapo Hyvärinen and Erkki Oja",
        "Independent Component Analysis: Algorithms and Applications",
        "2000, Neural Networks 13(4–5), 411–430",
        "https://doi.org/10.1016/S0893-6080(00)00026-5",
        "Statistical assumptions and identifiability in blind source separation.",
    ),
    "mmd": (
        "Arthur Gretton et al.",
        "A Kernel Two-Sample Test",
        "2012, JMLR 13, 723–773",
        "https://www.jmlr.org/papers/v13/gretton12a.html",
        "Mean embeddings and maximum mean discrepancy.",
    ),
    "keops": (
        "Benjamin Charlier, Jean Feydy, Joan Alexis Glaunès, François-David Collin, and Ghislain Durif",
        "Kernel Operations on the GPU, with Autodiff, without Memory Overflows",
        "2021, JMLR 22(74), 1–6",
        "https://www.jmlr.org/papers/v22/20-275.html",
        "Scalable reductions of large pairwise interaction matrices.",
    ),
}
GROUPS = {
    "image": ["mallat", "peyre", "skimage", "scipy", "strang"],
    "wavelets": ["mallat", "daubechies", "peyre", "donoho", "scipy"],
    "coding": ["cover", "shannon", "huffman", "mallat", "peyre"],
    "denoising": ["nlmreview", "donoho", "stein", "rof", "mallat", "peyre"],
    "nlm": ["nlm", "nlmreview", "mallat", "skimage", "stein"],
    "block": ["cai", "donoho", "stein", "mallat", "daubechies"],
    "optim": ["boyd", "prox", "fista", "chambollepock", "condat"],
    "inverse": ["rof", "prox", "fista", "chambollepock", "mallat", "boyd"],
    "completion": ["completion", "svt", "prox", "boyd", "strang"],
    "sparse": ["cs", "lasso", "prox", "fista", "mallat"],
    "dictionary": ["ksvd", "mairal", "mallat", "lasso", "prox"],
    "ml": ["esl", "strang", "sklearn", "boyd", "lasso"],
    "sgd": ["sgd", "adam", "esl", "boyd", "deep"],
    "deep": ["deep", "backprop", "adam", "sgd", "esl"],
    "ot": ["ot", "sinkhorn", "scaling", "benamou", "boyd"],
    "unbalanced": ["unbalanced", "ot", "sinkhorn", "benamou", "prox"],
    "semidiscrete": ["semidiscrete", "ot", "sinkhorn", "boyd", "sgd"],
    "shapes": ["canny", "snakes", "levelset", "chanvese", "peyre"],
    "geodesic": ["fastmarch", "dijkstra", "levelset", "isomap", "peyre"],
    "geometry": ["geometry", "floater", "desbrun", "strang", "peyre"],
    "fluids": ["stam", "bridson", "levelset", "scipy", "peyre"],
    "gaussian": ["galerne", "portilla", "mallat", "strang", "peyre"],
    "audio": ["bss", "ica", "mallat", "peyre", "scipy"],
    "diffusion": ["ddpm", "scoresde", "vincent", "hyvarinen", "deep", "adam"],
    "conformal": ["conformal", "vovk", "shafer", "lei", "cqr"],
    "graphical": ["graphlasso", "banerjee", "meinshausen", "lasso", "boyd", "esl"],
    "particles": ["mmd", "keops", "ot", "sgd", "deep"],
    "neuraltexture": ["gatys", "portilla", "galerne", "deep", "backprop"],
}
# Each entry supplies a topic-specific introduction and searchable vocabulary.
ROWS = """
introduction_3_image|Introduction to Image Processing|Foundations|Explore how arrays become images, how convolution removes fine detail, and how derivatives reveal boundaries. Compare spatial operations with their Fourier-domain counterparts to connect the formulas with visible changes in an image.|pixels convolution blur gradients FFT Fourier frequency filtering|image
introduction_4_fourier_wavelets|Fourier and Wavelet Approximation|Foundations|Represent the same image in Fourier and wavelet coordinates, then reconstruct it from a limited number of coefficients. The comparison reveals why localized edges favor multiscale representations while smooth oscillations favor Fourier modes.|FFT Fourier Haar wavelets approximation sparsity coefficients thresholding|wavelets
introduction_4_fourier_wavelets_interactive|Explore Fourier and Wavelets Interactively|Foundations|Use interactive controls to vary how many transform coefficients an image keeps. Compare reconstruction quality across Fourier and wavelet representations, and connect abrupt visual changes with the coefficients being discarded.|widgets sliders interactive FFT compression image wavelet approximation|wavelets
introduction_5_image_approx|Wavelet and DCT Image Compression|Foundations|Compare wavelet and discrete cosine representations for image approximation. Reconstruct images at different coefficient budgets to see how localization, block structure, and thresholding influence compression artifacts.|DCT cosine JPEG wavelets compression threshold coefficients artifacts|wavelets
introduction_6_elementary_fr|Le traitement numérique des images|Foundations|Découvrez comment représenter une image par une matrice, modifier son contraste et construire des filtres simples. Les expériences relient chaque opération numérique à son effet visuel et préparent l’étude des représentations multi-échelles.|français débutant image contraste filtrage bruit ondelettes matrices|image
wavelet_4_daubechies2d|Two-Dimensional Daubechies Wavelets|Wavelets & compression|Build a separable orthogonal wavelet transform from its analysis filters and reconstruct the image with the corresponding synthesis filters. Follow coefficients across scales to understand vanishing moments, localization, and perfect reconstruction.|Daubechies orthogonal wavelets filter bank multiresolution vanishing moments reconstruction|wavelets
multidim_2_volumetric|Volumetric Wavelet Processing|Wavelets & compression|Extend multiscale image processing to a three-dimensional volume. Examine slices and isosurfaces, then study how coefficient thresholding changes both the interior signal and the geometry of its level sets.|3D volume Haar wavelets isosurface medical imaging thresholding denoising|wavelets
coding_1_approximation|Image Approximation in Orthogonal Bases|Wavelets & compression|Compare linear and nonlinear approximations in several orthogonal bases. By keeping a fixed coefficient budget, we can measure how efficiently each basis represents smooth regions, textures, and edges.|orthogonal basis Fourier wavelet DCT compression approximation sparsity error|wavelets
coding_2_entropic|Entropy Coding and Compression|Wavelets & compression|Construct Huffman codes from symbol probabilities and check that encoding followed by decoding recovers the data exactly. Group symbols into blocks to see how average code lengths approach the entropy bound.|Huffman entropy Shannon information theory prefix code lossless compression source coding|coding
denoisingsimp_2b_linear_image|Linear Image Denoising|Image restoration|Study denoising as a tradeoff between suppressing random fluctuations and retaining useful detail. Compare linear filters in the Fourier domain and relate their frequency responses to the reconstructed images.|Wiener Gaussian blur linear filter Fourier noise restoration bias variance|denoising
denoisingwav_2_wavelet_2d|Wavelet Image Denoising|Image restoration|Remove noise by shrinking wavelet coefficients and reconstructing the image. Compare threshold choices to see how sparse representations separate meaningful structure from small noise-dominated coefficients.|wavelet shrinkage soft hard threshold Gaussian noise SNR PSNR denoising|denoising
denoisingwav_4_block|Wavelet Block Thresholding|Image restoration|Process neighboring wavelet coefficients together instead of thresholding each coefficient independently. The experiments show how shared local energy changes the balance between texture preservation and noise removal.|block threshold wavelets James Stein shrinkage local energy denoising|block
denoisingadv_6_nl_means|Nonlocal Means Denoising|Image restoration|Find similar image patches and average their center pixels using similarity weights. Explore the effects of patch size and search windows, then examine how a lower-dimensional patch representation changes the estimate.|nonlocal means NLM patch similarity PCA Gaussian weights texture denoising|nlm
denoisingadv_7_rankfilters|Rank Filters and Mathematical Morphology|Image restoration|Replace local pixel averages with order statistics such as medians and quantiles. Compare their behavior on noisy images and connect extreme ranks with erosion, dilation, and other morphological operations.|median quantile rank erosion dilation morphology salt pepper robust filtering|denoising
denoisingadv_9_sure|Stein Unbiased Risk Estimation|Image restoration|Choose denoising parameters by estimating reconstruction risk from the noisy observation itself. Compare the estimated risk with the error available in a controlled experiment and study where Stein’s identity enters the calculation.|SURE Stein unbiased risk Monte Carlo divergence threshold parameter selection denoising|denoising
inverse_2_deconvolution_variational|Variational Image Deconvolution|Inverse problems|Recover an image from a blurred observation by balancing data fidelity with a regularity penalty. Follow the optimization iterates and compare reconstructions to understand why direct inversion amplifies noise.|deconvolution inverse problem total variation Tikhonov Fourier regularization gradient|inverse
inverse_5_inpainting_sparsity|Sparse Image Inpainting|Inverse problems|Reconstruct missing pixels by promoting sparsity in a wavelet representation. Alternate consistency with the observed pixels and coefficient shrinkage, then inspect how the sampling pattern affects recovery.|inpainting missing pixels sparse wavelets ISTA regularization image reconstruction|inverse
inverse_6_trace_norm_dr|Matrix Completion with the Nuclear Norm|Inverse problems|Recover a low-rank matrix from a subset of its entries. Singular-value thresholding and Douglas–Rachford splitting provide a concrete way to combine a rank surrogate with exact agreement on observed values.|matrix completion low rank nuclear trace norm SVD Douglas Rachford missing data|completion
sparsity_6_l1_recovery|Sparse Recovery with ℓ¹ Minimization|Inverse problems|Investigate when a sparse vector can be recovered from fewer measurements than unknowns. Compare recovery outcomes across sparsity and sampling levels, and relate the experiments to the geometry of ℓ¹ minimization.|compressed sensing l1 basis pursuit sparse recovery phase transition undersampling|sparse
sparsity_4_dictionary_learning|Learning Sparse Dictionaries|Inverse problems|Learn an overcomplete collection of image atoms by alternating sparse coding and dictionary updates. Inspect the learned atoms and the objective history to connect the optimization problem with recurring local image structures.|dictionary learning sparse coding patches overcomplete atoms matrix factorization K SVD|dictionary
optim_1_gradient_descent|Gradient Descent|Optimization|Follow gradient descent on smooth objectives and apply it to image restoration. Changing the step size reveals the relationship between curvature, stability, and the speed at which the objective decreases.|gradient descent step size Lipschitz convex quadratic total variation optimization|optim
optim_1bis_gradient_descent|Gradient and Projected Gradient Methods|Optimization|Develop gradient descent for quadratic objectives, image denoising, and constrained inpainting. Explicit gradients and projections make it possible to check the operator adjoints and observe how step-size choices influence convergence.|Laurent Condat gradient projected descent Tikhonov inpainting quadratic constraints|optim
optim_2_condat_fb|Forward–Backward Splitting|Optimization|Recover a sparse signal by combining a smooth data-fidelity gradient with soft thresholding. Compare basic, relaxed, and accelerated iterations to see how different updates solve the same nonsmooth problem.|Laurent Condat proximal gradient ISTA FISTA forward backward sparse deconvolution|optim
optim_3_condat_fb_duality|Forward–Backward Splitting on the Dual|Optimization|Solve total-variation denoising through its dual formulation. Projection onto pixelwise Euclidean balls replaces a difficult primal proximity operator, while primal and dual energies provide a useful convergence diagnostic.|Laurent Condat dual total variation TV denoising FISTA primal dual gap projection|optim
optim_4_condat_dr|Douglas–Rachford Splitting|Optimization|Recover a sparse signal under exact linear constraints by alternating two proximity operators. Track feasibility and sparsity as the Douglas–Rachford iterates evolve, then repeat the experiment with a less sparse signal.|Laurent Condat Douglas Rachford proximal splitting basis pursuit sparse recovery constraints|optim
optim_5_condat_primal_dual|Chambolle–Pock Primal–Dual Splitting|Optimization|Implement a primal–dual algorithm for total-variation inpainting and denoising. The complete examples make the image gradient, its adjoint, and each proximity operator explicit, and check the result using feasibility and a primal–dual gap.|Chambolle Pock Laurent Condat primal dual TV inpainting denoising proximal gap|optim
optim_6_interior_points|Interior-Point Methods|Optimization|Approach constrained optimization from the interior of the feasible set using a logarithmic barrier. Follow Newton updates along the central path and inspect how the barrier parameter changes feasibility and optimality.|interior point Newton barrier central path linear programming convex optimization|optim
optim_7_noncvx_pro|Smooth Bilevel Sparse Regularization|Optimization|Reformulate a sparse regression objective through a smooth factorization and an inner linear solve. Check the reduced gradient numerically and compare the resulting optimization trajectory with direct sparse regression methods.|bilevel variable projection nonconvex sparse lasso factorization gradient L BFGS|optim
optimaltransp_1_linprog|Optimal Transport by Linear Programming|Optimal transport|Compute a minimum-cost coupling between discrete probability distributions. Inspect the transport plan and its marginals, and compare the geometry of the coupling with the cost that defines the optimization problem.|Kantorovich Monge Wasserstein linear programming coupling assignment transport plan|ot
optimaltransp_5_entropic|Entropic Optimal Transport|Optimal transport|Replace a linear transport problem with a smooth entropically regularized objective. Derive alternating scaling updates and explore how the regularization level changes the coupling, convergence, and transported data.|Sinkhorn entropy regularization coupling Wasserstein matrix scaling transport|ot
optimaltransp_5_entropic_full|Entropic Transport and Applications|Optimal transport|Use Sinkhorn scaling to compare and transform discrete distributions. Move from transport plans to applications, and examine how entropy changes the sharpness of a coupling as well as the numerical behavior of the iterations.|Sinkhorn optimal transport color transfer entropy barycenter coupling Wasserstein|ot
optimaltransp_6_entropic_adv|Advanced Sinkhorn Algorithms|Optimal transport|Explore refinements of entropic transport beyond basic alternating scaling. Inspect the effects of numerical stabilization and acceleration, and use marginal residuals to assess whether the computed coupling solves the intended problem.|Sinkhorn stabilization log domain acceleration entropy marginal residual optimal transport|ot
optimaltransp_7_semidiscrete|Semi-Discrete Optimal Transport|Optimal transport|Transport a continuous distribution toward a finite set of weighted points. Study the dual variables and the regions they induce, and connect geometric partitions with iterative optimization of the transport objective.|semi discrete Laguerre power diagram Voronoi stochastic gradient transport weights|semidiscrete
optimaltransp_8_unbalanced|Unbalanced Optimal Transport|Optimal transport|Allow transport to create or remove mass when exact marginal matching is inappropriate. Compare generalized scaling updates with the balanced case and inspect how marginal penalties change the resulting coupling.|unbalanced transport KL Kullback Leibler mass creation destruction entropy Sinkhorn|unbalanced
segmentation_1_edge_detection|Image Edge Detection|Geometry & segmentation|Locate image boundaries using derivatives at several scales. Compare gradient magnitudes, directional suppression, and level-set displays to distinguish a strong response from a precisely localized contour.|Canny edges gradient nonmaximum suppression Gaussian scale space segmentation|shapes
segmentation_2_snakes_param|Parametric Active Contours|Geometry & segmentation|Evolve a closed curve under smoothing and image-dependent forces. Track how curvature and external forces cooperate to attract the curve toward object boundaries, and inspect the role of initialization and step size.|snakes parametric active contour curvature curve evolution segmentation image forces|shapes
segmentation_3_snakes_levelset|Level-Set Active Contours|Geometry & segmentation|Represent a moving contour as the zero level of a function and evolve that function on a grid. Compare geometric and region-based forces, and examine how redistancing helps maintain a usable level-set representation.|level set snakes Chan Vese curvature redistancing segmentation PDE front propagation|shapes
fastmarching_0_implementing|Dijkstra and Fast Marching|Geometry & segmentation|Propagate a distance front through a weighted domain, first with a graph update and then with an eikonal update. Extract paths from the resulting distance map to see how local travel costs determine global geodesics.|Dijkstra fast marching eikonal shortest path geodesic distance front propagation|geodesic
shapes_7_isomap|Manifold Learning with Isomap|Geometry & segmentation|Recover a low-dimensional representation of a curved point cloud using neighborhood-graph distances. Compare geodesic and Euclidean distances, then use a spectral embedding to reveal the manifold’s intrinsic geometry.|Isomap manifold learning geodesic graph shortest path MDS dimensionality reduction swiss roll|geodesic
meshdeform_1_parameterization|Mesh Parameterization|Geometry & segmentation|Map a triangulated surface into the plane by fixing its boundary and solving for interior coordinates. Examine the resulting distortion and apply a texture to connect the linear system with a geometric visualization.|mesh parameterization harmonic barycentric Laplacian texture mapping boundary surface|geometry
meshproc_3_denoising|Mesh Denoising|Geometry & segmentation|Smooth noisy vertex positions using discrete differential operators on a triangulated surface. Compare the geometric effect of the updates and inspect how the choice of weights influences detail preservation and shrinkage.|mesh denoising Laplace Beltrami normals surface smoothing cotangent geometry|geometry
audio_2_separation|Audio Source Separation with Sparsity|Signals & graphics|Separate several audio sources from a smaller number of mixtures using sparse time-frequency structure. Examine spectrograms, estimate mixing directions, and listen to reconstructed signals to understand both the promise and limitations of masking.|audio sound blind source separation STFT spectrogram sparsity mixture microphone masking|audio
graphics_1_synthesis_gaussian|Gaussian Texture Synthesis|Signals & graphics|Generate random textures with prescribed second-order statistics using Fourier-domain operations. Compare samples with a reference texture to see what covariance captures and which structured features a Gaussian model cannot reproduce.|Gaussian texture covariance stationary random field Fourier phase synthesis spectrum|gaussian
graphics_5_fluids|Fluid Dynamics|Signals & graphics|Advect a velocity field and an image through a two-dimensional domain. Combine interpolation with an incompressibility projection, and inspect how these numerical ingredients create coherent fluid motion.|fluid Navier Stokes incompressible advection divergence projection semi Lagrangian simulation|fluids
ml_1_pca_nn|PCA, Nearest Neighbors, and Clustering|Machine learning|Explore a labeled dataset through principal components, neighborhood classifiers, and clustering. Compare low-dimensional visualizations with prediction results to separate the roles of representation, distance, and supervision.|PCA SVD principal components nearest neighbors kNN kmeans clustering iris classification|ml
ml_2_regression|Linear Regression, Ridge, and Lasso|Machine learning|Fit and compare linear regression models with quadratic and sparse penalties. Keep training and test data separate while inspecting coefficient paths and prediction errors, so that improved fit is not confused with improved generalization.|regression least squares ridge lasso prostate regularization cross validation SGD prediction|ml
ml_3_classification|Logistic and Kernel Classification|Machine learning|Turn linear scores into class probabilities and optimize a logistic loss. Move from linear to kernel decision boundaries, then inspect how regularization and representation affect the separation of the classes.|logistic regression classification kernel Gaussian decision boundary cross entropy sigmoid|ml
ml_4_sgd|Stochastic Gradient Descent|Machine learning|Replace full gradients with random sample updates and compare the resulting trajectories. Study step-size schedules and variance reduction to understand why noisy iterations can be effective for large learning problems.|SGD stochastic gradient variance reduction SAG learning rate logistic regression|sgd
ml_5_mlp|Multilayer Perceptrons|Machine learning|Build a small neural network and differentiate its loss through successive layers. Inspect the learned decision regions and the optimization history to connect nonlinear feature transformations with classification.|MLP multilayer perceptron neural network backpropagation activation classification|deep
ml_5bis_mlp_approx|Neural Networks for Function Approximation|Machine learning|Train a multilayer network to approximate a nonlinear function. Compare architectures and inspect where the fitted function differs from the target, relating expressivity to the optimization problem used to learn the weights.|neural network approximation PyTorch MLP activation width depth regression|deep
ml_6_nn|Deep Neural Networks and Automatic Differentiation|Machine learning|Implement forward and backward propagation, then compare manual derivatives with automatic differentiation. Vary the network architecture to see how depth and nonlinear activations influence both the fitted boundary and training dynamics.|deep neural network backpropagation autograd PyTorch gradient architecture classification|deep
ml_7_trend_filtering|Simplex Regression with Trend Filtering|Machine learning|Estimate a structured signal under simplex constraints and a trend-filtering penalty. Compare piecewise-polynomial behavior with unconstrained fitting, and inspect how the regularization parameter changes the recovered structure.|trend filtering simplex regression total variation piecewise polynomial convex CVXPY|inverse
ml_8_deep_texture_synthesis|Deep Texture Synthesis|Machine learning|Synthesize a texture by matching Gram matrices of pretrained convolutional features. Optimize the image pixels while keeping the feature extractor fixed, and compare this nonlinear representation with classical texture statistics.|texture synthesis VGG CNN Gram matrix style transfer PyTorch convolution deep features|neuraltexture
ml_9_graphical_lasso|Graphical Lasso|Machine learning|Infer sparse conditional dependencies by estimating a Gaussian precision matrix. Start with a known synthetic graph, follow a regularization path, and select a penalty using held-out likelihood rather than access to the true graph.|graphical lasso precision inverse covariance Gaussian graph conditional independence sparse|graphical
ml_10_particle_system|Interacting Particle Systems|Machine learning|Evolve a cloud of particles under pairwise kernel interactions and compare analytic and automatic gradients. Explore periodic boundaries and chunked computations to connect the continuum interpretation with a practical discrete simulation.|particles kernel MMD interaction gradient flow periodic boundary PyTorch chunking KeOps|particles
ml_11_conformal_prediction|Conformal Prediction|Machine learning|Construct prediction regions by comparing a candidate observation with observed conformity scores. Visualize the resulting sets and connect their coverage interpretation with exchangeability and the finite-sample ranking argument.|conformal prediction uncertainty coverage calibration exchangeability p value regression confidence|conformal
ml_12_diffusion_models|Denoising Diffusion Models|Machine learning|Study forward noising and reverse sampling on a one-dimensional Gaussian mixture. Because the exact score is available, we can compare numerical sampling with the target density before training a neural approximation by denoising score matching.|diffusion generative model score matching Langevin SDE Ornstein Uhlenbeck PyTorch denoising|diffusion
ml_12_diffusion_model_jax|Diffusion Models with JAX|Machine learning|Implement Gaussian-mixture diffusion and score learning with JAX and Flax. Compare exact-score sampling with a learned score field, using automatic differentiation and explicit random keys to make the stochastic computation reproducible.|JAX Flax Optax diffusion score matching Langevin generative SDE automatic differentiation|diffusion
"""
VARIANTS = {
    "ml_2_regression-BAD": "ml_2_regression",
    "ml_2_regression-SVG": "ml_2_regression",
    "ml_6_nn_OLD": "ml_6_nn",
    "optimaltransp_1_linprog-OLD": "optimaltransp_1_linprog",
}


def catalog():
    data = {}
    for row in ROWS.strip().splitlines():
        slug, title, category, intro, keywords, group = row.split("|")
        data[slug] = dict(
            slug=slug,
            title=title,
            category=category,
            intro=intro,
            keywords=keywords.split(),
            references=GROUPS[group],
        )
    for slug in data:
        if slug.startswith(("ml_5", "ml_6", "ml_8", "ml_12")):
            data[slug]["keywords"] += [
                "neural networks",
                "deep learning",
                "automatic differentiation",
            ]
        if slug.startswith("denoising"):
            data[slug]["keywords"] += [
                "image restoration",
                "noise removal",
                "filtering",
            ]
        if slug.startswith("optimaltransp"):
            data[slug]["keywords"] += [
                "optimal transportation",
                "probability distributions",
                "Wasserstein distance",
            ]
    for slug, canonical in VARIANTS.items():
        data[slug] = dict(data[canonical], slug=slug, canonical=canonical)
    return data

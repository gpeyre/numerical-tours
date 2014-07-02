<?
begin_toc();
toc_entry('Introduction', 'introduction');
toc_entry('Wavelet Processing', 'wavelet');
toc_entry('Approximation, Coding and Compression', 'coding');
toc_entry('Simple Denoising Methods', 'denoisingsimp');
toc_entry('Wavelet Denoising', 'denoisingwav');
toc_entry('Advanced Denoising Methods', 'denoisingadv');
toc_entry('Audio Processing', 'audio');
toc_entry('Higher Dimensional Signal Processing', 'multidim');
toc_entry('Computer Graphics', 'graphics');
toc_entry('Optimal Transport', 'optimaltransp');
toc_entry('Numerical Analysis', 'numerics');
toc_entry('Optimization', 'optim');
toc_entry('Variational Image Processing', 'variational');
toc_entry('Sparsity and Redundant Representations', 'sparsity');
toc_entry('Inverse Problems', 'inverse');
toc_entry('Compressive Sensing', 'cs');
toc_entry('Geodesic Processing', 'fastmarching');
toc_entry('Shapes', 'shapes');
toc_entry('Mesh Processing', 'meshproc');
toc_entry('Mesh Parameterization and Deformation', 'meshdeform');
toc_entry('Multiscale Mesh Processing', 'meshwav');
end_toc();

begin_tours('Introduction', 'introduction');
tour('introduction_1_basics', 'Basic Matlab/Scilab Instructions');
tour('introduction_2_signal', 'Introduction to Signal Processing');
tour('introduction_3_image', 'Introduction to Image Processing');
tour('introduction_4_fourier_wavelets', 'Image Approximation with Fourier and Wavelets');
tour('introduction_5_wavelets_2d', 'Image Processing with Wavelets');
tour('introduction_6_elementary_fr', 'Le traitement numérique des images');
end_tours();

begin_tours('Wavelet Processing', 'wavelet');
tour('wavelet_1_haar1d', '1-D Haar Wavelets');
tour('wavelet_2_haar2d', '2-D Haar Wavelets');
tour('wavelet_3_daubechies1d', '1-D Daubechies Wavelets');
tour('wavelet_4_daubechies2d', '2-D Daubechies Wavelets');
end_tours();

begin_tours('Approximation, Coding and Compression', 'coding');
tour('coding_1_approximation', 'Image Approximation with Orthogonal Bases');
tour('coding_2_entropic', 'Entropic Coding and Compression');
tour('coding_3_natural_images', 'Natural Images Statistics');
tour('coding_4_wavelet_compression', 'Image Compression with Wavelets');
tour('coding_5_watermarking', 'Wavelet Domain Image Watermarking');
end_tours();

begin_tours('Simple Denoising Methods', 'denoisingsimp');
tour('denoisingsimp_1_noise_models', 'Signal and Image Noise Models');
tour('denoisingsimp_2_linear', 'Image Denoising with Linear Methods');
end_tours();

begin_tours('Wavelet Denoising', 'denoisingwav');
tour('denoisingwav_1_wavelet_1d', 'Signal Denoising with Wavelets');
tour('denoisingwav_2_wavelet_2d', 'Image Denoising with Wavelets');
tour('denoisingwav_3_advanced', 'Advanced Wavelet Thresholdings');
tour('denoisingwav_4_block', 'Wavelet Block Thresholding');
tour('denoisingwav_5_data_dependent', 'Data Dependent Noise Models');
tour('denoisingwav_6_curvelets', 'Curvelet Denoising');
end_tours();

begin_tours('Advanced Denoising Methods', 'denoisingadv');
tour('denoisingadv_1_denoiseflow', 'Denoising by Sobolev and Total Variation Flows');
tour('denoisingadv_2_denoiseregul', 'Denoising by Sobolev and Total Variation Regularization');
tour('denoisingadv_3_chambollealgo', 'Total Variation Regularization with Chambolle Algorihtm');
tour('denoisingadv_4_median', 'Outliers and Median Denoiser');
tour('denoisingadv_5_mathmorph', 'Mathematical Morphology');
tour('denoisingadv_6_nl_means', 'Non Local Means');
tour('denoisingadv_7_rankfilters', 'Rank Filters for Image Processing');
tour('denoisingadv_8_bilateral', 'Bilateral Filtering');
tour('denoisingadv_9_sure', 'Stein Unbiased Risk Estimator');
end_tours();

begin_tours('Audio Processing', 'audio');
tour('audio_1_processing', 'Sound Processing with Short Time Fourier Transform');
tour('audio_2_separation', 'Source Separation with Sparsity');
end_tours();

begin_tours('Higher Dimensional Signal Processing', 'multidim');
tour('multidim_1_color', 'Color Image Processing');
tour('multidim_2_volumetric', 'Volumetric wavelet Data Processing');
tour('multidim_3_multispectral', 'Multi-spectral Imaging');
tour('multidim_4_video', 'Video Coding');
tour('multidim_5_opticalflow', 'Optical Flow Computation');
tour('multidim_6_tomography', 'Volumetric Radon Inversion');
tour('multidim_7_median', 'Color Image Denoising with Median Filtering');
end_tours();

begin_tours('Computer Graphics', 'graphics');
tour('graphics_1_synthesis_gaussian', 'Gaussian Models for Texture Synthesis');
tour('graphics_2_synthesis_wavelets', 'Texture Synthesis Using Wavelets');
tour('graphics_3_synthesis_diffusion', 'Texture Synthesis with PDEs');
tour('graphics_4_multiplicative_cascade', 'Multiplicative Cascade Synthesis of Signals and Textures');
tour('graphics_5_fluids', 'Fluid Dynamics');
tour('graphics_6_patches', 'Texture Synthesis and Inpainting using Patch Projections');
tour('graphics_7_shape_shading', 'Shape From Shading');
tour('graphics_8_dyntextures', 'Stationary Dynamic Texture Synthesis');
end_tours();

begin_tours('Optimal Transport', 'optimaltransp');
tour('optimaltransp_1_linprog', 'Optimal Transport with Linear Programming');
tour('optimaltransp_2_benamou_brenier', 'Optimal Transport with Benamou-Brenier Algorithm');
tour('optimaltransp_3_matching_1d', 'Optimal Transport in 1-D');
tour('optimaltransp_4_matching_sliced', 'Sliced Optimal Transport');
end_tours();

begin_tours('Numerical Analysis', 'numerics');
tour('numerics_1_wavelet_compression', 'Wavelet Compression of Integral Operators');
end_tours();

begin_tours('Optimization', 'optim');
tour('optim_1_gradient_descent', 'Gradient Descent Methods');
tour('optim_2_newton', 'Newton Method');
tour('optim_3_cgs', 'Conjugate Gradient');
end_tours();

begin_tours('Variational Image Processing', 'variational');
tour('variational_1_heat', 'Heat Diffusion');
tour('variational_2_edge_detection', 'Edge Detection');
tour('variational_3_image_separation', 'Cartoon+Texture Variational Image Decomposition');
tour('variational_4_snakes_param', 'Active Contours using Parameteric Curves');
tour('variational_5_snakes_levelset', 'Active Contours using Level Sets');
tour('variational_6_structure_tensor', 'Anisotropic Diffusion with Structure Tensor');
tour('variational_7_convex_segmentation', 'Convex Region-Based Image Segmentation');
end_tours();

begin_tours('Sparsity and Redundant Representations', 'sparsity');
tour('sparsity_1_homotopy', 'L1 Minimization with Homotopy');
tour('sparsity_2_matching_pursuit', 'Sparse Spikes Deconvolution with Matching Pursuits');
tour('sparsity_3_gabor', 'Sparse Representation in a Gabor Dictionary');
tour('sparsity_4_dictionary_learning', 'Dictionary Learning');
tour('sparsity_4bis_dictionary_learning_denoising', 'Dictionary Learning for Denoising');
tour('sparsity_5_sudoku', 'Sudoku using POCS and Sparsity');
end_tours();

begin_tours('Inverse Problems', 'inverse');
tour('inverse_1_sparse_spikes', 'Sparse 1D Deconvolution');
tour('inverse_2_deconvolution_variational', 'Image Deconvolution using Variational Method');
tour('inverse_3_deconvolution_sparsity', 'Image Deconvolution using Sparse Regularization');
tour('inverse_4_inpainting_variational', 'Inpainting using Variational Regularization');
tour('inverse_5_inpainting_sparsity', 'Inpainting using Sparse Regularization');
tour('inverse_6_primal_dual', 'Total Variation Regularization using Primal Dual Schemes');
tour('inverse_7_nl_inpainting', 'Inpainting with NL-means');
tour('inverse_8_tomography', 'Reconstruction from Partial Tomography Measurements');
tour('inverse_9_tomography_sobsparse', 'Tomography Inversion using Tikhonov and Sparse Regularization');
tour('inverse_9b_gfb', 'Generalized Forward-Backward Splitting');
end_tours();

begin_tours('Compressive Sensing', 'cs');
tour('cs_1_sparse_signal', 'Compressed Sensing of Sparse Signals');
tour('cs_2_fourier', 'Reconstruction from Compressive Fourier Measurements');
tour('cs_3_images', 'Compressed Sensing of Images');
tour('cs_4_matrix_completion', 'Matrix Completion with Nuclear Norm Minimization');
tour('cs_5_l1_recovery', 'Performance of Sparse Recovery Using L1 Minimization');
end_tours();

begin_tours('Geodesic Processing', 'fastmarching');
tour('fastmarching_0_implementing', 'Dijkstra and Fast Marching Algorithms');
tour('fastmarching_1_2d', 'Fast Marching in 2D');
tour('fastmarching_2_3d', 'Fast Marching in 3D');
tour('fastmarching_3_anisotropy', 'Anisotropic Fast Marching');
tour('fastmarching_4_mesh', 'Geodesic Mesh Processing');
tour('fastmarching_4bis_geodesic_mesh', 'Geodesic Distance Computation on 3-D Meshes');
tour('fastmarching_5_sampling_2d', 'Geodesic Farthest Point Sampling');
tour('fastmarching_6_sampling_surf', 'Geodesic Surface Remeshing');
tour('fastmarching_7_sampling_compr', 'Geodesic Triangulation for Image Compression');
tour('fastmarching_8_segmentation', 'Geodesic Segmentation');
tour('fastmarching_9_heuristics', 'Heuristically Driven Front Propagation');
end_tours();

begin_tours('Shapes', 'shapes');
tour('shapes_1_bendinginv_2d', 'Geodesic Bending Invariants for Shapes');
tour('shapes_2_bendinginv_3d', 'Geodesic Bending Invariants for Surfaces');
tour('shapes_3_bendinginv_landmarks', 'Geodesic Bending Invariants with Landmarks');
tour('shapes_4_shape_matching', 'Shape Correspondence with Fast Marching');
tour('shapes_5_geodesic_descriptors', 'Shape Retrieval with Geodesic Descriptors');
tour('shapes_6_medialaxis', 'Geodesic Medial Axsis');
tour('shapes_7_isomap', 'Manifold Learning with Isomap');
end_tours();

begin_tours('Mesh Processing', 'meshproc');
tour('meshproc_1_basics_2d', 'Basics About 2D Triangulation');
tour('meshproc_2_basics_3d', 'Basics About 3D Meshes');
tour('meshproc_3_denoising', 'Mesh Denoising');
tour('meshproc_4_fourier', 'Fourier on Meshes');
tour('meshproc_5_pde', 'PDEs on Meshes');
tour('meshproc_6_volumetric', 'Volumetric Meshes');
end_tours();

begin_tours('Mesh Parameterization and Deformation', 'meshdeform');
tour('meshdeform_1_parameterization', 'Mesh Parameterization');
tour('meshdeform_2_parameterization_sphere', 'Spherical Mesh Parameterization');
tour('meshdeform_3_flattening', 'Spectral Mesh Flattening');
tour('meshdeform_4_barycentric', 'Generalized Barycentric Coordinates for Warpping');
tour('meshdeform_5_deformation', 'Mesh Deformation');
end_tours();

begin_tours('Multiscale Mesh Processing', 'meshwav');
tour('meshwav_1_subdivision_curves', 'Subdivision Curves');
tour('meshwav_2_subdivision_surfaces', 'Subdivision Surfaces');
tour('meshwav_3_simplification', 'Mesh Simplification');
tour('meshwav_4_haar_sphere', 'Spherical Haar Wavelets');
tour('meshwav_5_wavelets', 'Wavelet Transform on 3D Meshes');
end_tours();

?>

function perform_index_generation_php(filename)

if nargin<1
    filename = '../html/index_tours.php';
end

fid = fopen(filename, 'wt');
if fid<=0
    error('Unable to open file');
end

list_ext = {...
       {'introduction'  'Introduction'} ...
       {'wavelet' 'Wavelet Processing'} ... 
       {'coding'  'Approximation, Coding and Compression'} ... 
       {'denoisingsimp'  'Simple Denoising Methods'} ...  
       {'denoisingwav'  'Wavelet Denoising'} ...  
       {'denoisingadv'  'Advanced Denoising Methods'} ...  
       {'audio' 'Audio Processing'} ...  
       {'multidim' 'Higher Dimensional Signal Processing'} ... 
       {'graphics' 'Computer Graphics'} ... 
       {'numerics' 'Numerical Analysis'} ... 
       {'optim' 'Optimization'} ...
       {'variational' 'Variational Image Processing'}  ... 
       {'sparsity' 'Sparsity and Redundant Representations'} ... 
       {'inverse' 'Inverse Problems'} ... 
       {'cs' 'Compressive Sensing'} ... 
       {'fastmarching' 'Geodesic Processing'} ...
       {'shapes' 'Shapes'} ...
       {'meshproc' 'Mesh Processing'} ... 
       {'meshdeform' 'Mesh Parameterization and Deformation'} ... 
       {'meshwav' 'Multiscale Mesh Processing'} ... 
    };

pr = @(x)fprintf(fid,[x '\n']);
prL = @()fprintf(fid, '\n');

pr('<?');

%%% GENERATE TOC %%%
pr('begin_toc();');
for iext = 1:length(list_ext)
    ext = list_ext{iext}{1};
    tit = list_ext{iext}{2};    
    pr(['toc_entry(''' tit ''', ''' ext ''');']);
end
pr('end_toc();');
prL();

%%% GENERATE SECTIONS %%%
for iext = 1:length(list_ext)
    ext = list_ext{iext}{1};
    tit = list_ext{iext}{2};
    a = dir([ext '_*.m']);
    
    pr(['begin_tours(''' tit ''', '''  ext ''');']);
    for k=1:length(a)
        tourname = a(k).name;
        fidt = fopen(tourname);
        L = fgets(fidt);
        L = strtrim(strrep(L, '%% ', ''));
        fclose(fidt);
        pr(['tour(''' tourname(1:end-2) ''', '''  L  ''');']);
    end
	pr('end_tours();');
    prL();
end


pr('?>');

fclose(fid);



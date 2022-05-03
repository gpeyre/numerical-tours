#include "fm2dAniso.h"

// PBR - MKR 20220503
// old prototype in /Applications/MATLAB_R2021b.app/extern/include/matrix.h
// LIBMMWMATRIX_PUBLISHED_API_EXTERN_C mxArray *mxCreateNumericArray(mwSize ndim, const mwSize *dims, mxClassID classid, mxComplexity flag);
// type of param 2 not correct: old = const mwSize *dims   new = const int dims[3]
#undef mxCreateNumericArray
LIBMMWMATRIX_PUBLISHED_API_EXTERN_C mxArray * mxCreateNumericArray(mwSize ndim, const int dims[3], mxClassID classid, mxComplexity flag);
LIBMMWMATRIX_PUBLISHED_API_EXTERN_C mxArray * mxCreateNumericArray(mwSize ndim, const int dims[3], mxClassID classid, mxComplexity flag){
    return mxCreateNumericArray_730(ndim,(const mwSize *)dims,classid,flag);
}
 

// PBR - MKR 20220503
// old prototype in /Applications/MATLAB_R2021b.app/extern/include/matrix.h
// LIBMMWMATRIX_PUBLISHED_API_EXTERN_C const mwSize *mxGetDimensions(const mxArray *pa);
// type of return not correct: old = const mwSize *   new = const int*
#undef mxGetDimensions
LIBMMWMATRIX_PUBLISHED_API_EXTERN_C const int *mxGetDimensions(const mxArray *pa);
LIBMMWMATRIX_PUBLISHED_API_EXTERN_C const int *mxGetDimensions(const mxArray *pa) {
    return (const int*)mxGetDimensions_730(pa);
}

// PBR - MKR 20220503
// old prototype in /Applications/MATLAB_R2021b.app/extern/include/matrix.h
// LIBMMWMATRIX_PUBLISHED_API_EXTERN_C int mxSetDimensions(mxArray *pa, const mwSize *pdims, mwSize ndims);
// type of param 2 not correct: old = const mwSize *dims   new = const int *dims
#undef mxSetDimensions
LIBMMWMATRIX_PUBLISHED_API_EXTERN_C int mxSetDimensions(mxArray *pa, const int *pdims, mwSize ndims);
LIBMMWMATRIX_PUBLISHED_API_EXTERN_C int mxSetDimensions(mxArray *pa, const int *pdims, mwSize ndims) {
    return mxSetDimensions_730(pa, (const mwSize *)pdims, ndims);
}



void mexFunction(	int nlhs, mxArray *plhs[], 
					int nrhs, const mxArray*prhs[] ) 
{
	//------------------------------------------------------------------
	/* retrive arguments */
	//------------------------------------------------------------------
	if( (nrhs!=3) && (nrhs!=4) )
		mexErrMsgTxt("3 or 4 arguments are required.");
	if( mxGetNumberOfDimensions(prhs[1])!= 4 )
		mexErrMsgTxt("HERE T must be a 2D x 2x2 tensor field of symmetric definite matrices.");
	//------------------------------------------------------------------
	// First argument : spacing and dimensions
	const int* dim_h = mxGetDimensions(prhs[0]);
    if ( (dim_h[0]!=2) || (dim_h[1]!=1) )
	  mexErrMsgTxt("Library error: h must be a 2x1 array list.");
	hx = mxGetPr(prhs[0])[0]; hy = mxGetPr(prhs[0])[1];
	hx2 = hx*hx; hy2 = hy*hy;
    hxhy = hx*hy;
	hx2hy2 = hx*hx*hy*hy;
	hx2_plus_hy2 = hx*hx + hy*hy;
    sqrt_hx2_plus_hy2 = sqrt(hx2_plus_hy2);
    nx = mxGetDimensions(prhs[1])[0];
	ny = mxGetDimensions(prhs[1])[1];
    if( (mxGetDimensions(prhs[1])[2] != 2) || (mxGetDimensions(prhs[1])[3] != 2) )
        mexErrMsgTxt("T must be a 2D x 2x2 tensor field of symmetric definite matrices.");
	Nx = nx+2; Ny = ny+2;
	size = Nx*Ny; nxny = nx*ny;
	//------------------------------------------------------------------
	// Second argument : Anisotropy matrices 
	T = mxGetPr(prhs[1]);
	//------------------------------------------------------------------
	// Third argument : start points
	start_points = mxGetPr(prhs[2]);
	nb_start_points = mxGetN(prhs[2]);
	//------------------------------------------------------------------
    Dmax = 0.0;
	//------------------------------------------------------------------
    // forth argument : max length to reach 
    if(nrhs==4)
    	Dmax = (float) mxGetScalar(prhs[3]);
	//==================================================================
	// Outputs
	int dims[2] = {Nx,Ny};
	// First output : minimal action map
	plhs[0] = mxCreateNumericArray(2, dims, mxSINGLE_CLASS, mxREAL );
	U = (float*) mxGetPr(plhs[0]);
	// Second output : dUx
	plhs[1] = mxCreateNumericArray(2, dims, mxSINGLE_CLASS, mxREAL );
	dUx = (float*) mxGetPr(plhs[1]);
	// Third output : dUy
	plhs[2] = mxCreateNumericArray(2, dims, mxSINGLE_CLASS, mxREAL );
	dUy = (float*) mxGetPr(plhs[2]);
	// Fourth output : Voronoi diagramm
	plhs[3] = mxCreateNumericArray(2, dims, mxINT16_CLASS, mxREAL );
	V = (short*) mxGetPr(plhs[3]);
	// Fifth output : Euclidien Distance
	plhs[4] = mxCreateNumericArray(2, dims, mxSINGLE_CLASS, mxREAL );
	L = (float*) mxGetPr(plhs[4]);
	//==================================================================
	InitializeNeighborhoods();
	//------------------------------------------------------------------
	InitializeArrays();
	//------------------------------------------------------------------
	InitializeOpenHeap();
	//------------------------------------------------------------------
	RunPropagation();
    //==================================================================
    resize();
    dims[0] = Nx-2; dims[1] = Ny-2;
    mxSetDimensions(plhs[0], dims, 2);
    mxSetDimensions(plhs[1], dims, 2);
    mxSetDimensions(plhs[2], dims, 2);
    mxSetDimensions(plhs[3], dims, 2);
    mxSetDimensions(plhs[4], dims, 2);
    //==================================================================
	DELETEARRAY(NeighborhoodLarge);
    DELETEARRAY(Simplicies1);
    DELETEARRAY(Simplicies2);
    DELETEARRAY(signsX);
    DELETEARRAY(signsY);
    DELETEARRAY(h1);
    DELETEARRAY(h2);
    DELETEARRAY(h1_h2);
    DELETEARRAY(h22);
    DELETEARRAY(q_gradient);
    DELETEARRAY(S);
	DELETEARRAY(M1);DELETEARRAY(M2);DELETEARRAY(M3);
    DELETEARRAY(Trial);
    DELETEARRAY(Tree);
    return;
}

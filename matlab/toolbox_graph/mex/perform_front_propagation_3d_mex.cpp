/*=================================================================
% perform_front_propagation_3d - perform a Fast Marching front propagation.
%
%   OLD : [D,S] = perform_front_propagation_2d(W,start_points,end_points,nb_iter_max,H);
%	[D,S,Q] = perform_front_propagation_3d(W,start_points,end_points,nb_iter_max, H, L);
%
%   'D' is a 2D array containing the value of the distance function to seed.
%	'S' is a 2D array containing the state of each point : 
%		-1 : dead, distance have been computed.
%		 0 : open, distance is being computed but not set.
%		 1 : far, distance not already computed.
%	'W' is the weight matrix (inverse of the speed).
%	'start_points' is a 3 x num_start_points matrix where k is the number of starting points.
%	'H' is an heuristic (distance that remains to goal). This is a 2D matrix.
%   
%   Copyright (c) 2004 Gabriel Peyré
*=================================================================*/

#include "perform_front_propagation_3d.h"
#include "mex.h"


void mexFunction(	int nlhs, mxArray *plhs[], 
					int nrhs, const mxArray*prhs[] ) 
{ 
	/* retrive arguments */
	if( nrhs<4 ) 
		mexErrMsgTxt("4 - 7 input arguments are required."); 
	if( nlhs<1 ) 
		mexErrMsgTxt("1, 2 or 3 output arguments are required."); 

	// first argument : weight list
	if( mxGetNumberOfDimensions(prhs[0])!= 3 )
		mexErrMsgTxt("W must be a 3D array.");
	n = mxGetDimensions(prhs[0])[0];
	p = mxGetDimensions(prhs[0])[1];
	q = mxGetDimensions(prhs[0])[2];
	W = mxGetPr(prhs[0]);
	// second argument : start_points
	start_points = mxGetPr(prhs[1]);
	int tmp = mxGetM(prhs[1]); 
	nb_start_points = mxGetN(prhs[1]);
	if( nb_start_points==0 || tmp!=3 )
		mexErrMsgTxt("start_points must be of size 3 x nb_start_poins."); 
	// third argument : end_points
	end_points = mxGetPr(prhs[2]);
	tmp = mxGetM(prhs[2]); 
	nb_end_points = mxGetN(prhs[2]);
	if( nb_end_points!=0 && tmp!=3 )
		mexErrMsgTxt("end_points must be of size 3 x nb_end_poins."); 
	// argument 4 : nb_iter_max
	nb_iter_max = (int) *mxGetPr(prhs[3]);
	// argument 5 : heuristic
	if( nrhs>=5 )
	{
		H = mxGetPr(prhs[4]);
		if( mxGetM(prhs[4])==0 && mxGetN(prhs[4])==0 )
			H=NULL;
		if( H!=NULL && (mxGetDimensions(prhs[4])[0]!=n || mxGetDimensions(prhs[4])[1]!=p || mxGetDimensions(prhs[4])[2]!=q) )
			mexErrMsgTxt("H must be of size n x p x q."); 
	}
	else
		H = NULL;
	// argument 6 : constraint map
	if( nrhs>=6 )
	{
		L = mxGetPr(prhs[5]);
		if( L!=NULL && (mxGetDimensions(prhs[5])[0]!=n || mxGetDimensions(prhs[5])[1]!=p || mxGetDimensions(prhs[5])[2]!=q) )
			mexErrMsgTxt("L must be of size n x p x q."); 
	}
	else
		L = NULL;
	// argument 7: value list
	if( nrhs>=7 )
	{
		values = mxGetPr(prhs[6]);
		if( mxGetM(prhs[6])==0 && mxGetN(prhs[6])==0 )
			values=NULL;
		if( values!=NULL && (mxGetM(prhs[6])!=nb_start_points || mxGetN(prhs[6])!=1) )
			mexErrMsgTxt("values must be of size nb_start_points x 1."); 
	}
	else
		values = NULL;
		
	// first ouput : distance
	int dims[3] = {n,p,q};
	plhs[0] = mxCreateNumericArray(3, dims, mxDOUBLE_CLASS, mxREAL );
	D = mxGetPr(plhs[0]);
	// second output : state
	if( nlhs>=2 )
	{
		plhs[1] = mxCreateNumericArray(3, dims, mxDOUBLE_CLASS, mxREAL );
		S = mxGetPr(plhs[1]);
	}
	else
	{
		S = new double[n*p*q];
	}
	// third output : index
	if( nlhs>=3 )
	{
		plhs[2] = mxCreateNumericArray(3, dims, mxDOUBLE_CLASS, mxREAL );
		Q = mxGetPr(plhs[2]);
	}
	else
	{
		Q = new double[n*p*q];
	}

	
	// launch the propagation
	perform_front_propagation_3d();

	if( nlhs<2 )
		GW_DELETEARRAY(S);		
	if( nlhs<3 )
		GW_DELETEARRAY(Q);

	return;
}
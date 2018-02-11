/*=================================================================
% perform_front_propagation_2d - perform a Fast Marching front propagation.
%
%   [D,S,Q] = perform_front_propagation_2d(W,start_points,end_points,nb_iter_max,H,L,values);
%
%   'D' is a 2D array containing the value of the distance function to seed.
%	'S' is a 2D array containing the state of each point : 
%		-1 : dead, distance have been computed.
%		 0 : open, distance is being computed but not set.
%		 1 : far, distance not already computed.
%	Q is the index of the closest point.
%	'W' is the weight matrix (inverse of the speed).
%	'start_points' is a 2 x num_start_points matrix where k is the number of starting points.
%	'H' is an heuristic (distance that remains to goal). This is a 2D matrix.
%	L is a constraint matrix, points will be considered only if their current distance is less than L.
%   
%   Copyright (c) 2004 Gabriel Peyré
*=================================================================*/

#include "perform_front_propagation_2d.h"
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
	n = mxGetM(prhs[0]); 
	p = mxGetN(prhs[0]);
	W = mxGetPr(prhs[0]);
	// second argument : start_points
	start_points = mxGetPr(prhs[1]);
	int tmp = mxGetM(prhs[1]); 
	nb_start_points = mxGetN(prhs[1]);
	if( nb_start_points==0 || tmp!=2 )
		mexErrMsgTxt("start_points must be of size 2 x nb_start_poins."); 
	// third argument : end_points
	end_points = mxGetPr(prhs[2]);
	tmp = mxGetM(prhs[2]); 
	nb_end_points = mxGetN(prhs[2]);
	if( nb_end_points!=0 && tmp!=2 )
		mexErrMsgTxt("end_points must be of size 2 x nb_end_poins."); 
	//  argument 4: nb_iter_max
	nb_iter_max = (int) *mxGetPr(prhs[3]);
	//  argument 5: heuristic
	if( nrhs>=5 )
	{
		H = mxGetPr(prhs[4]);
		if( mxGetM(prhs[4])==0 && mxGetN(prhs[4])==0 )
			H=NULL;
		if( H!=NULL && (mxGetM(prhs[4])!=n || mxGetN(prhs[4])!=p) )
			mexErrMsgTxt("H must be of size n x p.");  
	}
	else
		H = NULL;
	// argument 6: constraint map
	if( nrhs>=6 )
	{
		L = mxGetPr(prhs[5]);
		if( mxGetM(prhs[5])==0 && mxGetN(prhs[5])==0 )
			H=NULL;
		if( L!=NULL && (mxGetM(prhs[5])!=n || mxGetN(prhs[5])!=p) )
			mexErrMsgTxt("L must be of size n x p."); 
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
	plhs[0] = mxCreateDoubleMatrix(n, p, mxREAL); 
	D = mxGetPr(plhs[0]);
	// second output : state
	if( nlhs>=2 )
	{
		plhs[1] = mxCreateDoubleMatrix(n, p, mxREAL); 
		S = mxGetPr(plhs[1]);
	}
	else
	{
		S = new double[n*p];
	}
	// third output : index
	if( nlhs>=3 )
	{
		plhs[2] = mxCreateDoubleMatrix(n, p, mxREAL); 
		Q = mxGetPr(plhs[2]);
	}
	else
	{
		Q = new double[n*p];
	} 
	
	// launch the propagation
	perform_front_propagation_2d(); 
	

	if( nlhs<2 )
		GW_DELETEARRAY(S);		
	if( nlhs<3 )
		GW_DELETEARRAY(Q);
	return;
}

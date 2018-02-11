/*=================================================================
% perform_circular_front_propagation_2d - perform a Fast Marching front propagation.
%
%   [D,S] = perform_circular_front_propagation_2d(W,start_points,end_points,center_point,nb_iter_max,H);
%
%   'D' is a 2D array containing the value of the distance function to seed.
%	'S' is a 2D array containing the state of each point : 
%		-1 : dead, distance have been computed.
%		 0 : open, distance is being computed but not set.
%		 1 : far, distance not already computed.
%	'W' is the weight matrix (inverse of the speed).
%	'start_points' is a 2 x num_start_points matrix where k is the number of starting points.
%	'H' is an heuristic (distance that remains to goal). This is a 2D matrix.
%   
%   Copyright (c) 2004 Gabriel Peyré
*=================================================================*/

#include "mex.h"
#include "perform_front_propagation_2d.h"

double* center_point = NULL;

bool callback_intert_node(int i, int j, int ii, int jj)
{
	// point must be on the right of the center
	if( i>=center_point[0] && ii>=(int)center_point[0] )
	{
		if( j>=center_point[1] && jj<(int)center_point[1] )
			return false;
		if( j<center_point[1] && jj>=(int)center_point[1] )
			return false;
	}
	return true;
}

void mexFunction(	int nlhs, mxArray *plhs[], 
				 int nrhs, const mxArray*prhs[] ) 
{ 
	/* retrive arguments */
	if( nrhs<5 ) 
		mexErrMsgTxt("5 or 6 input arguments are required."); 
	if( nlhs<1 ) 
		mexErrMsgTxt("1 or 2 output arguments are required."); 

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
	// fourth argument : center_point
	center_point = mxGetPr(prhs[3]);
	// fifth argument : nb_iter_max
	nb_iter_max = (int) *mxGetPr(prhs[4]);
	// sixth argument : heuristic
	if( nrhs==6 )
	{
		H = mxGetPr(prhs[5]);
		if( mxGetM(prhs[5])!=n || mxGetN(prhs[5])!=p )
			mexErrMsgTxt("H must be of size n x p."); 
	}
	else
		H = NULL;

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

	// launch the propagation
	perform_front_propagation_2d(callback_intert_node);

	if( nlhs<2 )
		GW_DELETEARRAY(S);
	return;
}

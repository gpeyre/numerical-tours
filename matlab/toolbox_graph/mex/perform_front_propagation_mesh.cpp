/*=================================================================
% perform_front_propagation_mesh - perform a Fast Marching front propagation on a 3D mesh.
%
%   [D,S,Q] = perform_front_propagation_mesh(vertex, faces, W,start_points,end_points, nb_iter_max,H,L, values, dmax);
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

#include <math.h>
#include "config.h"
#include <algorithm>
#include <map>
#include <vector>
#include <list>
#include <string>
#include <iostream>
#include <fstream>
#include <string.h>
using std::string;
using std::cerr;
using std::cout;
using std::endl;

#include "mex.h"
#include "gw/gw_core/GW_Config.h"
#include "gw/gw_core/GW_MathsWrapper.h"
#include "gw/gw_geodesic/GW_GeodesicMesh.h"
using namespace GW;


inline void display_message(const char* mess, int v)
{
	char str[128];
	sprintf(str, mess, v);
	mexWarnMsgTxt(str);
}



double* vertex = NULL;
int nverts = -1; 
double* faces = NULL;
int nfaces = -1; 
double* start_points = NULL;
int nstart = -1;
double* end_points = NULL;
int nend = -1;
double* H = NULL;	// heuristic
double* L = NULL;	// bound on current distance
double* Ww = NULL;	// weight
int niter_max = -1;
double dmax = 1e9;
double* values = NULL;
// outputs 
double* D = NULL;	// distance
double* S = NULL;	// state
double* Q = NULL;	// nearest neighbor

#define faces_(k,i) faces[k+3*i]
#define vertex_(k,i) vertex[k+3*i]


GW_Float WeightCallback(GW_GeodesicVertex& Vert)
{
	GW_U32 i = Vert.GetID();
	return Ww[i];
}

GW_Bool StopMarchingCallback( GW_GeodesicVertex& Vert )
{
	// check if the end point has been reached
	GW_U32 i = Vert.GetID();
//	display_message("ind %d",i );
//	display_message("dist %f",Vert.GetDistance() );
	if( Vert.GetDistance()>dmax )
		return true;
	for( int k=0; k<nend; ++k )
		if( end_points[k]==i )
			return true;
	return false;
}
int nbr_iter = 0;
GW_Bool InsersionCallback( GW_GeodesicVertex& Vert, GW_Float rNewDist )
{
	// check if the distance of the new point is less than the given distance
	GW_U32 i = Vert.GetID();
	bool doinsersion = nbr_iter<=niter_max;
	if( L!=NULL )
		doinsersion = doinsersion && (rNewDist<L[i]);
	nbr_iter++;
	return doinsersion;
}
GW_Float HeuristicCallback( GW_GeodesicVertex& Vert )
{
	// return the heuristic distance
	GW_U32 i = Vert.GetID();
	return H[i];
}


void mexFunction(	int nlhs, mxArray *plhs[], 
				 int nrhs, const mxArray*prhs[] ) 
{ 
	nbr_iter = 0;
	/* retrive arguments */
	if( nrhs<6 ) 
		mexErrMsgTxt("6 or 7 input arguments are required."); 
	if( nlhs<1 ) 
		mexErrMsgTxt("1 or 2 output arguments are required."); 

	// arg1 : vertex
	vertex = mxGetPr(prhs[0]);
	nverts = mxGetN(prhs[0]); 
	if( mxGetM(prhs[0])!=3 )
		mexErrMsgTxt("vertex must be of size 3 x nverts."); 
	// arg2 : faces
	faces = mxGetPr(prhs[1]);
	nfaces = mxGetN(prhs[1]);
	if( mxGetM(prhs[1])!=3 )
		mexErrMsgTxt("face must be of size 3 x nfaces."); 
	// arg3 : W
	Ww = mxGetPr(prhs[2]);
	int m = mxGetM(prhs[2]);
	if( m!=nverts )
		mexErrMsgTxt("W must be of same size as vertex."); 
	// arg4 : start_points
	start_points = mxGetPr(prhs[3]);
	nstart = mxGetM(prhs[3]);
	// arg5 : end_points
	end_points = mxGetPr(prhs[4]);
	nend = mxGetM(prhs[4]);
	// arg6 : niter_max
	niter_max = (int) *mxGetPr(prhs[5]);
	// arg7 : H
	if( nrhs>=7 )
	{
		H = mxGetPr(prhs[6]);
		int m =mxGetM(prhs[6]);
		if( m>0 && m!=nverts )
			mexErrMsgTxt("H must be of size nverts."); 
		if( m==0 )
			H = NULL;
	}
	else
	{
		H = NULL;
	}
	// arg8 : L
	if( nrhs>=8 )
	{
		L = mxGetPr(prhs[7]);
		int m =mxGetM(prhs[7]);
		if( m>0 && mxGetM(prhs[7])!=nverts )
			mexErrMsgTxt("L must be of size nverts."); 
		if( m==0 )
			L = NULL;
	}
	else
		L = NULL;
		
	// argument 9: value list
	if( nrhs>=9 )
	{
		values = mxGetPr(prhs[8]);
		if( mxGetM(prhs[8])==0 && mxGetN(prhs[8])==0 )
			values=NULL;
		if( values!=NULL && (mxGetM(prhs[8])!=nstart || mxGetN(prhs[8])!=1) )
			mexErrMsgTxt("values must be of size nb_start_points x 1."); 
	}
	else
		values = NULL;
	// argument 10: dmax
	if( nrhs>=9 )
		dmax = *mxGetPr(prhs[9]);
	else
		dmax = 1e9;


	// first ouput : distance
	plhs[0] = mxCreateDoubleMatrix(nverts, 1, mxREAL); 
	D = mxGetPr(plhs[0]);
	// second output : state
	plhs[1] = mxCreateDoubleMatrix(nverts, 1, mxREAL); 
	S = mxGetPr(plhs[1]);
	// second output : segmentation
	plhs[2] = mxCreateDoubleMatrix(nverts, 1, mxREAL); 
	Q = mxGetPr(plhs[2]);

	// create the mesh
	GW_GeodesicMesh Mesh;
	Mesh.SetNbrVertex(nverts);
	for( int i=0; i<nverts; ++i )
	{
		GW_GeodesicVertex& vert = (GW_GeodesicVertex&) Mesh.CreateNewVertex();
		vert.SetPosition( GW_Vector3D(vertex_(0,i),vertex_(1,i),vertex_(2,i)) );
		Mesh.SetVertex(i, &vert);
	}
	Mesh.SetNbrFace(nfaces);
	for( int i=0; i<nfaces; ++i )
	{
		GW_GeodesicFace& face = (GW_GeodesicFace&) Mesh.CreateNewFace();
		GW_Vertex* v1 = Mesh.GetVertex((int) faces_(0,i)); GW_ASSERT( v1!=NULL );
		GW_Vertex* v2 = Mesh.GetVertex((int) faces_(1,i)); GW_ASSERT( v2!=NULL );
		GW_Vertex* v3 = Mesh.GetVertex((int) faces_(2,i)); GW_ASSERT( v3!=NULL );
		face.SetVertex( *v1,*v2,*v3 );
		Mesh.SetFace(i, &face);
	}
	Mesh.BuildConnectivity();

	// set up fast marching	
	Mesh.ResetGeodesicMesh();
	for( int i=0; i<nstart; ++i )
	{
		GW_GeodesicVertex* v = (GW_GeodesicVertex*) Mesh.GetVertex((GW_U32) start_points[i]);
		GW_ASSERT( v!=NULL );
		Mesh.AddStartVertex( *v );
	}
	Mesh.SetUpFastMarching();
	Mesh.RegisterWeightCallbackFunction( WeightCallback );
	Mesh.RegisterForceStopCallbackFunction( StopMarchingCallback );
	Mesh.RegisterVertexInsersionCallbackFunction( InsersionCallback );
	if( H!=NULL )
		Mesh.RegisterHeuristicToGoalCallbackFunction( HeuristicCallback );
	// initialize the distance of the starting points
	if( values!=NULL )
	for( int i=0; i<nstart; ++i )
	{
		GW_GeodesicVertex* v = (GW_GeodesicVertex*) Mesh.GetVertex((GW_U32) start_points[i]);
		GW_ASSERT( v!=NULL );
		v->SetDistance( values[i] );
	}
	
	// perform fast marching
//	display_message("itermax=%d", niter_max);
	Mesh.PerformFastMarching();

	// output result
	for( int i=0; i<nverts; ++i )
	{
		GW_GeodesicVertex* v = (GW_GeodesicVertex*) Mesh.GetVertex((GW_U32) i);
		GW_ASSERT( v!=NULL );
		D[i] = v->GetDistance();
		S[i] = v->GetState();
		GW_GeodesicVertex* v1 = v->GetFront();
		if( v1==NULL )
			Q[i] = -1;
		else
			Q[i] = v1->GetID();
	}
	

	return;
}

/*=================================================================
% perform_front_propagation_2d - perform a Fast Marching front propagation.
%
%   [D,S,Q] = perform_front_propagation_2d(W,start_points,end_points,nb_iter_max, H, L, values);
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
%   Copyright (c) 2004 Gabriel Peyr�
*=================================================================*/

/* 
Warning : there is a bug in the fibheap data.
Quand on procède à un remplacement de données 
(ligne 302-304), cela induit une erreur de tri dans la structure de l'arbre vu que le 
code de fib. ne copie pas les données mais elles sont passées par pointeur. Autrement dit, 
l'appel de la fonction fh_replacedata() ne change rien dans ta version actuelle. 
Il faut absolument l'appeler avec un autre point de donnée. 
*/

// select to test or not to test (debug purpose)
// #define CHECK_HEAP check_heap(i,j,k);
#ifndef CHECK_HEAP
	#define CHECK_HEAP
#endif
// error display
// #define ERROR_MSG(a) mexErrMsgTxt(a)
#ifndef ERROR_MSG
	#define ERROR_MSG(a) 
#endif
// #define WARN_MSG(a)  mexWarnMsgTxt(a) 
#ifndef WARN_MSG
	#define WARN_MSG(a)
#endif

#include "perform_front_propagation_2d.h"
#include "fheap/fib.h"
#include "fheap/fibpriv.h"

#define kDead -1
#define kOpen 0
#define kFar 1

/* Global variables */
int n;			// width
int p;			// height
double* D = NULL;
double* S = NULL;
double* W = NULL;
double* Q = NULL;
double* start_points = NULL;
double* end_points = NULL;
double* values = NULL;
double* H = NULL;
double* L = NULL;
int nb_iter_max = 100000;
int nb_start_points = 0;
int nb_end_points = 0;
fibheap_el** heap_pool = NULL;

#define ACCESS_ARRAY(a,i,j) a[(i)+n*(j)]
#define D_(i,j) ACCESS_ARRAY(D,i,j)
#define S_(i,j) ACCESS_ARRAY(S,i,j)
#define W_(i,j) ACCESS_ARRAY(W,i,j)
#define H_(i,j) ACCESS_ARRAY(H,i,j)
#define Q_(i,j) ACCESS_ARRAY(Q,i,j)
#define L_(i,j) ACCESS_ARRAY(L,i,j)
#define heap_pool_(i,j) ACCESS_ARRAY(heap_pool,i,j)
#define start_points_(i,k) start_points[(i)+2*(k)]
#define end_points_(i,k) end_points[(i)+2*(k)]

struct point
{
	point( int ii, int jj )
	{ i = ii; j = jj; }
	int i,j;
};
typedef std::vector<point*> point_list;

inline 
bool end_points_reached(const int i, const int j )
{
	for( int k=0; k<nb_end_points; ++k )
	{
		if( i==((int)end_points_(0,k)) && j==((int)end_points_(1,k)) )
			return true;
	}
	return false;
}

inline 
int compare_points(void *x, void *y)
{
	point& a = *( (point*) x );
	point& b = *( (point*) y );
	if( H==NULL )
		return cmp( D_(a.i,a.j), D_(b.i,b.j) );
	else
		return cmp( D_(a.i,a.j)+H_(a.i,a.j), D_(b.i,b.j)+H_(b.i,b.j) );
}


// test the heap validity
void check_heap( int i, int j )
{
	for( int x=0; x<n; ++x )
		for( int y=0; y<p; ++y )
		{
			if( heap_pool_(x,y)!=NULL )
			{
				point& pt = * (point*)heap_pool_(x,y)->fhe_data;
				if( H==NULL )
				{
					if( D_(i,j)>D_(pt.i,pt.j) )
						ERROR_MSG("Problem with heap.\n");
				}
				else
				{
					if( D_(i,j)+H_(i,j)>D_(pt.i,pt.j)+H_(pt.i,pt.j) )
						ERROR_MSG("Problem with heap.\n");
				}
			}
		}
}



void perform_front_propagation_2d(T_callback_intert_node callback_insert_node)
{
	// create the Fibonacci heap
	struct fibheap* open_heap = fh_makeheap();
	fh_setcmp(open_heap, compare_points);

	double h = 1.0/n;
	
	// initialize points
	for( int i=0; i<n; ++i )
	for( int j=0; j<p; ++j )
	{
		D_(i,j) = GW_INFINITE;
		S_(i,j) = kFar;
		Q_(i,j) = -1;
	}

	// record all the points
	heap_pool = new fibheap_el*[n*p]; 
	memset( heap_pool, NULL, n*p*sizeof(fibheap_el*) );

	// inialize open list
	point_list existing_points;
	for( int k=0; k<nb_start_points; ++k )
	{
		int i = (int) start_points_(0,k);
		int j = (int) start_points_(1,k);

		if( D_( i,j )==0 )
			ERROR_MSG("start_points should not contain duplicates.");

		point* pt = new point( i,j );
		existing_points.push_back( pt );			// for deleting at the end
		heap_pool_(i,j) = fh_insert( open_heap, pt );			// add to heap
		if( values==NULL ) 
			D_( i,j ) = 0;
		else
			D_( i,j ) = values[k];
		S_( i,j ) = kOpen;
		Q_(i,j) = k;
	}

	// perform the front propagation
	int num_iter = 0;
	bool stop_iteration = GW_False;
	while( !fh_isempty(open_heap) && num_iter<nb_iter_max && !stop_iteration )
	{
		num_iter++;

		// current point
		point& cur_point = * ((point*) fh_extractmin( open_heap ));
		int i = cur_point.i;
		int j = cur_point.j;
		heap_pool_(i,j) = NULL;
		S_(i,j) = kDead;
		stop_iteration = end_points_reached(i,j);
		
		/*
		char msg[200];
		sprintf(msg, "Cool %f", Q_(i,j) );
		WARN_MSG( msg ); 
		*/
		
		CHECK_HEAP;

		// recurse on each neighbor
		int nei_i[4] = {i+1,i,i-1,i};
		int nei_j[4] = {j,j+1,j,j-1};
		for( int k=0; k<4; ++k )
		{
			int ii = nei_i[k];
			int jj = nei_j[k];
			bool bInsert = true;
			if( callback_insert_node!=NULL )
				bInsert = callback_insert_node(i,j,ii,jj);
			// check that the contraint distance map is ok
			if( ii>=0 && jj>=0 && ii<n && jj<p && bInsert )
			{
				double P = h/W_(ii,jj);
				// compute its neighboring values
				double a1 = GW_INFINITE;
				int k1 = -1;
				if( ii<n-1 )
				{
					bool bParticipate = true;
					if( callback_insert_node!=NULL )
						bParticipate = callback_insert_node(ii,jj,ii+1,jj);
					if( bParticipate )
					{
						a1 = D_(ii+1,jj);
						k1 = Q_(ii+1,jj);
					}
				}
				if( ii>0 )
				{
					bool bParticipate = true;
					if( callback_insert_node!=NULL )
						bParticipate = callback_insert_node(ii,jj,ii-1,jj);
					if( bParticipate )
					{
						if( D_(ii-1,jj)<a1 )
							k1 = Q_(ii-1,jj);
						a1 = GW_MIN( a1, D_(ii-1,jj) );
					}
				}
				double a2 = GW_INFINITE;
				int k2 = -1;
				if( jj<p-1 )
				{

					bool bParticipate = true;
					if( callback_insert_node!=NULL )
						bParticipate = callback_insert_node(ii,jj,ii,jj+1);
					if( bParticipate )
					{
						a2 = D_(ii,jj+1);
						k2 = Q_(ii,jj+1);
					}
				}
				if( jj>0 )
				{
					bool bParticipate = true;
					if( callback_insert_node!=NULL )
						bParticipate = callback_insert_node(ii,jj,ii,jj-1);
					if( bParticipate )
					{
						if( D_(ii,jj-1)<a2 )
							k2 = Q_(ii,jj-1);
						a2 = GW_MIN( a2, D_(ii,jj-1) );
					}
				}
				if( a1>a2 )	// swap so that a1<a2
				{
					double tmp = a1; a1 = a2; a2 = tmp;
					int tmpi = k1; k1 = k2; k2 = tmpi;
				}
				// update its distance
				// now the equation is   (a-a1)^2+(a-a2)^2 = P, with a >= a2 >= a1.
				double A1 = 0;
				if( 2*P*P > (a2-a1)*(a2-a1) ) // TO CHECK
				{
					double delta = 2*P*P-(a2-a1)*(a2-a1);
					A1 = (a1+a2+sqrt(delta))/2.0;
				}
				else
					A1 = a1 + P;
				if( ((int) S_(ii,jj)) == kDead )
				{
					// check if action has change. Should not happen for FM
					// if( A1<D_(ii,jj) )
					//	WARN_MSG("The update is not monotone");
#if 1
					if( A1<D_(ii,jj) )	// should not happen for FM
					{
						D_(ii,jj) = A1;
						// update the value of the closest starting point
						//if( GW_ABS(a1-A1)<GW_ABS(a2-A1) && k1>=0  )
							Q_(ii,jj) = k1;
						//else
						//	Q_(ii,jj) = k2;
						//Q_(ii,jj) = Q_(i,j);
					}
#endif
				}
				else if( ((int) S_(ii,jj)) == kOpen )
				{
					// check if action has change.
					if( A1<D_(ii,jj) )
					{
						D_(ii,jj) = A1;
						// update the value of the closest starting point
						//if( GW_ABS(a1-A1)<GW_ABS(a2-A1) && k1>=0  )
							Q_(ii,jj) = k1;
						//else
						//	Q_(ii,jj) = k2;
						//Q_(ii,jj) = Q_(i,j);
						// Modify the value in the heap
						fibheap_el* cur_el = heap_pool_(ii,jj);
						if( cur_el!=NULL )
							fh_replacedata( open_heap, cur_el, cur_el->fhe_data );	// use same data for update
						else
							ERROR_MSG("Error in heap pool allocation."); 
					}
				}
				else if( ((int) S_(ii,jj)) == kFar )
				{
					if( D_(ii,jj)!=GW_INFINITE )
						ERROR_MSG("Distance must be initialized to Inf");
					if( L==NULL || A1<=L_(ii,jj) )
					{
						S_(ii,jj) = kOpen;
						// distance must have change.
						D_(ii,jj) = A1;
						// update the value of the closest starting point
						//if( GW_ABS(a1-A1)<GW_ABS(a2-A1) && k1>=0 )
							Q_(ii,jj) = k1;
						//else
						//	Q_(ii,jj) = k2;
						//Q_(ii,jj) = Q_(i,j);
						// add to open list
						point* pt = new point(ii,jj);
						existing_points.push_back( pt );
						heap_pool_(ii,jj) = fh_insert( open_heap, pt );			// add to heap	
					}
				}
				else 
					ERROR_MSG("Unkwnown state."); 
					
			}	// end switch
		}		// end for
	}			// end while

//				char msg[200];
//				sprintf(msg, "Cool %f", Q_(100,100) );
//				 WARN_MSG( msg ); 

	// free heap
	fh_deleteheap(open_heap);
	// free point pool
	for( point_list::iterator it = existing_points.begin(); it!=existing_points.end(); ++it )
		GW_DELETE( *it );
	// free fibheap pool
	GW_DELETEARRAY(heap_pool);
}
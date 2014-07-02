
/*------------------------------------------------------------------------------*/
/** 
 *  \file   GW_Toolkit.h
 *  \brief  Definition of class \c GW_Toolkit
 *  \author Gabriel Peyré
 *  \date   4-27-2003
 */ 
/*------------------------------------------------------------------------------*/

#ifndef _GW_TOOLKIT_H_
#define _GW_TOOLKIT_H_

#include "../gw_core/GW_Config.h"
#include "../gw_geodesic/GW_GeodesicMesh.h"
#include "../gw_geodesic/GW_GeodesicPath.h"
#include "../gw_geodesic/GW_VoronoiMesh.h"
#include "../gw_core/GW_Mesh.h"
#include "../gw_toolkit/GW_ASELoader.h"
#include "../gw_toolkit/GW_PLYLoader.h"
#include "../gw_toolkit/GW_VRMLLoader.h"
#include "../gw_toolkit/GW_OFFLoader.h"
#include "../gw_toolkit/GW_OBJLoader.h"
#include "../gw_toolkit/GW_GeodesicDisplayer.h"
#include "../gw_geodesic/GW_GeometryAtlas.h"
#include "GW_OpenGLHelper.h"

namespace GW {

/*------------------------------------------------------------------------------*/
/** 
 *  \class  GW_Toolkit
 *  \brief  Generic helper functions.
 *  \author Gabriel Peyré
 *  \date   4-27-2003
 *
 *  Some usefull functions to manipulate a mesh and it's base parametrization.
 */ 
/*------------------------------------------------------------------------------*/

class GW_Toolkit
{

public:

	#define GW_TOOLKIT_VECTOR_SCALING 0.3

	GW_Toolkit();
	static GW_I32 LoadFileList( const char* file, T_StringList& FileList );
	GW_I32 LoadFile( const char* file_name ); 
	void SetUpMarching( GW_GeodesicVertex* pStartVertex = NULL );
	void ComputePath( GW_Bool bRandomizeVertex = GW_True );

	void AddFurthestPoint( GW_U32 nNbr = 1, GW_Bool bUseRandomStartVertex = GW_False,
						   GW_Bool bAssignColor = GW_True, GW_U32 nPrintEach = 1 );

	GW_GeodesicMesh& GetGeodesicMesh() { return Mesh_; }
	GW_VoronoiMesh& GetVoronoiMesh() { return VoronoiMesh_; }
	GW_GeodesicDisplayer& GetGeodesicDisplayer() { return Displayer_; }
	GW_BasicDisplayer& GetCoarseDisplayer() { return CoarseDisplayer_; }

	GW_GeodesicPath& GetGeodesicPath() { return GeodesicPath_; }

	void SetNbrStartVertex( GW_U32 nNbrStartVertex ) { nNbrStartVertex_=nNbrStartVertex; }
	GW_U32 AddStartVertex() { nNbrStartVertex_++; return nNbrStartVertex_; }
	GW_U32 RemoveStartVertex() { if( nNbrStartVertex_>0 ) nNbrStartVertex_--; return nNbrStartVertex_; }

	GW_GeodesicVertex* GetEndVertex()
	{ return pEndVertex_; }
	static T_GeodesicVertexList& GetStartVertex()
	{ return StartVertexList_; }

	static GW_Vector3D GetRandomColor();

	/** geometry cells helpers */
	static void Display( GW_GeometryCell& gc, GW_Float boundary_thick = 0, GW_Float half_boundary_thick = 0, 
						GW_Vector3D boundary_color = GW_Vector3D(0,0,0), GW_Vector3D half_boundary_color = GW_Vector3D(1,0,0), GW_Vector3D fill_color = GW_Vector3D(0.5f,0.5f,0.5f) );
	static void Display( GW_GeometryAtlas& ga, GW_Float boundary_thick = 0, GW_Float half_boundary_thick = 0, 
						GW_Vector3D boundary_color = GW_Vector3D(0,0,0), GW_Vector3D half_boundary_color = GW_Vector3D(1,0,0), GW_Vector3D fill_color = GW_Vector3D(0.5f,0.5f,0.5f) );

private:

	GW_GeodesicMesh Mesh_;
	GW_VoronoiMesh VoronoiMesh_;

	GW_GeodesicDisplayer Displayer_;
	GW_BasicDisplayer CoarseDisplayer_;

	GW_GeodesicPath GeodesicPath_;

	GW_GeodesicVertex* pEndVertex_;
	static T_GeodesicVertexList StartVertexList_;

	GW_U32 nNbrStartVertex_;


	#define GW_TOOLKIT_NBR_COLOR 18
	static const GW_Float rRandomColor_[GW_TOOLKIT_NBR_COLOR][3];

};

} // End namespace GW


#endif // _GW_TOOLKIT_H_


///////////////////////////////////////////////////////////////////////////////
//  Copyright (c) Gabriel Peyré
///////////////////////////////////////////////////////////////////////////////
//                               END OF FILE                                 //
///////////////////////////////////////////////////////////////////////////////

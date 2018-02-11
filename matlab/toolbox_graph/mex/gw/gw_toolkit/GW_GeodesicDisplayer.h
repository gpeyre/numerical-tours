
/*------------------------------------------------------------------------------*/
/** 
 *  \file   GW_GeodesicDisplayer.h
 *  \brief  Definition of class \c GW_GeodesicDisplayer
 *  \author Gabriel Peyré
 *  \date   4-10-2003
 */ 
/*------------------------------------------------------------------------------*/

#ifndef _GW_GEODESICDISPLAYER_H_
#define _GW_GEODESICDISPLAYER_H_

#include "../gw_core/GW_Config.h"
#include "../gw_geodesic/GW_GeodesicMesh.h"
#include "../gw_geodesic/GW_GeodesicPath.h"
#include "GW_BasicDisplayer.h"
#include "../gw_geodesic/GW_VoronoiMesh.h"

namespace GW {

/*------------------------------------------------------------------------------*/
/** \name a map of GW_Vector3D */
/*------------------------------------------------------------------------------*/
//@{
typedef std::map<GW_U32, GW_Vector3D> T_VertexColorMap;
typedef T_VertexColorMap::iterator IT_VertexColorMap;
typedef T_VertexColorMap::reverse_iterator RIT_VertexColorMap;
typedef T_VertexColorMap::const_iterator CIT_VertexColorMap;
typedef T_VertexColorMap::const_reverse_iterator CRIT_VertexColorMap;
//@}


/*------------------------------------------------------------------------------*/
/** 
 *  \class  GW_GeodesicDisplayer
 *  \brief  Display information about geodesic computations.
 *  \author Gabriel Peyré
 *  \date   4-10-2003
 *
 *  Use geodesic value to draw the mesh.
 */ 
/*------------------------------------------------------------------------------*/

class GW_GeodesicDisplayer: public GW_BasicDisplayer
{

public:

    /*------------------------------------------------------------------------------*/
    /** \name Constructor and destructor */
    /*------------------------------------------------------------------------------*/
    //@{
    GW_GeodesicDisplayer();
    //@}

	virtual void ComputeColor( GW_Vertex& pVert, float* color );
	virtual void SetUpDraw(GW_Mesh& Mesh);

	void DisplayPath( GW_GeodesicPath& CurPath, GW_Vector3D& Color = GW_Vector3D(1,1,1), GW_Float rLineWidth = 2.0f );
	void DisplayVoronoiPoints( GW_VoronoiMesh& MeshBuilder, GW_Vector3D& Color = GW_Vector3D(0,0,0), GW_Float rPointSize = 6.0f );
	void DisplayGeodesicBoundaries( GW_VoronoiMesh& VoronoiMesh, GW_Vector3D& Color = GW_Vector3D(0,0,0), GW_Float rLineWidth = 2.0f );

	void AddFrontColor( GW_Vertex& VertFront, GW_Vector3D& Color );
	void RemoveFrontColor( GW_Vertex& VertFront );
	void ResetFrontColor();

	GW_Vector3D GetStreamColor( GW_Float val, GW_Float max );
	void SetColorStreamFixedRadius(	GW_Float val );
	GW_Float GetMaxGeodesicDistance()
	{ return aMaxValue_[kGeodesicDistance]; }
	GW_Float GetColorStreamFixedRadius()
	{ return rColorStreamFixedRadius_; }

private:

	T_VertexColorMap VertexColorMap_;


	/** color strem management */
	#define COLOR_STREAM_SIZE 100
	GW_Vector3D ColorStream_[COLOR_STREAM_SIZE+1];
	GW_Float rColorStreamFixedRadius_;

};

} // End namespace GW

#ifdef GW_USE_INLINE
    #include "GW_GeodesicDisplayer.inl"
#endif


#endif // _GW_GEODESICDISPLAYER_H_


///////////////////////////////////////////////////////////////////////////////
//  Copyright (c) Gabriel Peyré
///////////////////////////////////////////////////////////////////////////////
//                               END OF FILE                                 //
///////////////////////////////////////////////////////////////////////////////

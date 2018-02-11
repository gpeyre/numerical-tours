/*------------------------------------------------------------------------------*/
/** 
 *  \file   GW_GeodesicDisplayer.cpp
 *  \brief  Definition of class \c GW_GeodesicDisplayer
 *  \author Gabriel Peyré
 *  \date   4-10-2003
 */ 
/*------------------------------------------------------------------------------*/


#ifdef GW_SCCSID
    static const char* sccsid = "@(#) GW_GeodesicDisplayer.cpp(c) Gabriel Peyré2003";
#endif // GW_SCCSID

#include "stdafx.h"
#include "GW_GeodesicDisplayer.h"
#include <GL/gl.h>

#ifndef GW_USE_INLINE
    #include "GW_GeodesicDisplayer.inl"
#endif

#include "GW_Toolkit.h"

using namespace GW;

/*------------------------------------------------------------------------------*/
// Name : GW_GeodesicDisplayer constructor
/**
 *  \author Gabriel Peyré
 *  \date   4-10-2003
 * 
 *  Constructor.
 */
/*------------------------------------------------------------------------------*/
GW_GeodesicDisplayer::GW_GeodesicDisplayer()
:	rColorStreamFixedRadius_ (-1)
{
	/* init the color stream */
	for( GW_U32	i=0; i<COLOR_STREAM_SIZE+1; ++i )
		ColorStream_[i] = GW_Toolkit::GetRandomColor();
}

/*------------------------------------------------------------------------------*/
// Name : GW_GeodesicDisplayer::GetStreamColor
/**
 *  \param  val [GW_Float] value (in [0,max])
 *  \param  max [GW_Float] max value
 *  \return [GW_Vector3D] Color from the flow
 *  \author Gabriel Peyré
 *  \date   9-29-2003
 * 
 *  Get a color from the stream, ouh yeah !
 */
/*------------------------------------------------------------------------------*/
GW_Vector3D GW_GeodesicDisplayer::GetStreamColor( GW_Float val, GW_Float max )
{
	val = val/max;
	GW_U32 nb = (GW_U32) floor(val*COLOR_STREAM_SIZE);
	GW_Float res = val - (nb*1.0f)/(COLOR_STREAM_SIZE*1.0f);
	GW_Vector3D& c1 = ColorStream_[(GW_U32) fmod( nb, COLOR_STREAM_SIZE)];
	GW_Vector3D& c2 = ColorStream_[(GW_U32) fmod( nb,COLOR_STREAM_SIZE)+1];
	res = res*COLOR_STREAM_SIZE;
	return c1*(1-res)+c2*res;
}


/*------------------------------------------------------------------------------*/
// Name : GW_GeodesicDisplayer::SetUpDraw
/**
 *  \param  Mesh [GW_Mesh&] The mesh.
 *  \author Gabriel Peyré
 *  \date   4-10-2003
 * 
 *  Set up the range for the geodesic distance.
 */
/*------------------------------------------------------------------------------*/
void GW_GeodesicDisplayer::SetUpDraw(GW_Mesh& Mesh)
{
	GW_BasicDisplayer::SetUpDraw( Mesh );
	GW_ASSERT( strcmp(Mesh.GetClassName().c_str(), "GW_GeodesicMesh")==0 );
	aMinValue_[kGeodesicDistance] = GW_INFINITE;
	aMaxValue_[kGeodesicDistance] = -GW_INFINITE;
	for( GW_U32 i = 0; i<Mesh.GetNbrVertex(); ++i )
	{
		GW_GeodesicVertex* pVert = (GW_GeodesicVertex*) Mesh.GetVertex(i);
		if( pVert->GetDistance()<GW_INFINITE*0.5 )
		{
			if( pVert->GetDistance()<aMinValue_[kGeodesicDistance] )
				aMinValue_[kGeodesicDistance] = pVert->GetDistance();
			if( pVert->GetDistance()>aMaxValue_[kGeodesicDistance] )
				aMaxValue_[kGeodesicDistance] = pVert->GetDistance();	
		}
	}
}

/*------------------------------------------------------------------------------*/
// Name : GW_GeodesicDisplayer::AddFrontColor
/**
 *  \param  VertFront [GW_Vertex&] The vertex from wich the front started.
 *  \param  Color [GW_Vector3D&] The color.
 *  \author Gabriel Peyré
 *  \date   4-10-2003
 * 
 *  Add a color to display a given front.
 */
/*------------------------------------------------------------------------------*/
void GW_GeodesicDisplayer::AddFrontColor( GW_Vertex& VertFront, GW_Vector3D& Color )
{
	GW_U32 nID = VertFront.GetID();
	VertexColorMap_[nID] = Color;
}

/*------------------------------------------------------------------------------*/
// Name : GW_GeodesicDisplayer::RemoveFrontColor
/**
 *  \param  VertFront [GW_Vertex&] The front.
 *  \author Gabriel Peyré
 *  \date   4-10-2003
 * 
 *  Remove a colorization.
 */
/*------------------------------------------------------------------------------*/
void GW_GeodesicDisplayer::RemoveFrontColor( GW_Vertex& VertFront )
{
	VertexColorMap_.erase( VertFront.GetID() );
}

/*------------------------------------------------------------------------------*/
// Name : GW_GeodesicDisplayer::ResetFrontColor
/**
 *  \author Gabriel Peyré
 *  \date   4-12-2003
 * 
 *  Remove all colors.
 */
/*------------------------------------------------------------------------------*/
void GW_GeodesicDisplayer::ResetFrontColor()
{
	VertexColorMap_.clear();
}


/*------------------------------------------------------------------------------*/
// Name : GW_GeodesicDisplayer::ComputeColor
/**
*  \param  pVert [GW_Vertex&] The vertex.
*  \param  color [GLfloat*] The color.
*  \author Gabriel Peyré
*  \date   4-3-2003
* 
*  Assign a color to a vertex according to it's curv info.
*/
/*------------------------------------------------------------------------------*/
void GW_GeodesicDisplayer::ComputeColor( GW_Vertex& Vert, float* color )
{
	GW_GeodesicVertex& GeoVert = (GW_GeodesicVertex&) Vert;
	GW_BasicDisplayer::ComputeColor( Vert, color );
	if( bProrieties[kGeodesicDistance] )
	{
		if( GeoVert.GetDistance()<GW_INFINITE/2 )
		{
			float dist = (float) GW_SCALE_01(GeoVert.GetDistance(), aMinValue_[kGeodesicDistance], aMaxValue_[kGeodesicDistance]);
			GW_ASSERT(GeoVert.GetFront()!=NULL);
			GW_U32 nID = GeoVert.GetFront()->GetID();
			if( VertexColorMap_.find(nID)!=VertexColorMap_.end() )
			{
				GW_Vector3D colvect = VertexColorMap_[nID];
				color[0] = (GLfloat) colvect[0]*(1-dist);
				color[1] = (GLfloat) colvect[1]*(1-dist);
				color[2] = (GLfloat) colvect[2]*(1-dist);
			}
			else
			{
				/* no colorisation was provided for this front */
				color[0] = color[1] = color[2] = 1-dist;
			}
			if( bProrieties[kGeodesicDistanceStreamColor] )
			{
				GW_Float max = aMaxValue_[kGeodesicDistance];
				if( rColorStreamFixedRadius_>0 )
					max = rColorStreamFixedRadius_;
				GW_Vector3D colvect = this->GetStreamColor( GeoVert.GetDistance(), max );
				color[0] = (GLfloat) colvect[0];
				color[1] = (GLfloat) colvect[1];
				color[2] = (GLfloat) colvect[2];
			}
		}
	}
	if( bProrieties[kMarchingState] )
	{

		if( GeoVert.GetState()==GW_GeodesicVertex::kAlive )
		{
			GW_Vector3D TheColor(1,0,0);	// default color = red
			GW_ASSERT(GeoVert.GetFront()!=NULL);
			GW_U32 nID = GeoVert.GetFront()->GetID();
			if( VertexColorMap_.find(nID)!=VertexColorMap_.end() )
				TheColor = VertexColorMap_[nID];
			color[0] = (float) TheColor[0]*0.5f;
			color[1] = (float) TheColor[1]*0.5f;
			color[2] = (float) TheColor[2]*0.5f;
		}
		if( GeoVert.GetState()==GW_GeodesicVertex::kDead )
		{
			GW_Vector3D TheColor(1,0,0);	// default color = red
			GW_ASSERT(GeoVert.GetFront()!=NULL);
			GW_U32 nID = GeoVert.GetFront()->GetID();
			if( VertexColorMap_.find(nID)!=VertexColorMap_.end() )
				TheColor = VertexColorMap_[nID];
			color[0] = (float) TheColor[0];
			color[1] = (float) TheColor[1];
			color[2] = (float) TheColor[2];
		}
	}
	if( bProrieties[kVertexParametrization] )
	{
#if 1
		GW_U32 nNbrParameter = 0;
		GW_Vector3D Color;
		for( GW_U32 i=0; i<3; ++i )
		{
			GW_Float rParam = 0;
			GW_VoronoiVertex* pParamVert = GeoVert.GetParameterVertex( i, rParam );
			if( pParamVert!=NULL )
			{
				GW_U32 nID = pParamVert->GetBaseVertex()->GetID();
				if( VertexColorMap_.find(nID)!=VertexColorMap_.end() )
				{
					Color += VertexColorMap_[nID]*rParam;
					nNbrParameter++;
				}
			}
		}
		color[0] = (float) Color[0];
		color[1] = (float) Color[1];
		color[2] = (float) Color[2];
		Color /= nNbrParameter;
#else
		for( GW_U32 i=0; i<3; ++i )
		{
			GW_Float rParam;
			GW_VoronoiVertex* pParamVert = GeoVert.GetParameterVertex( i, rParam );
			if( pParamVert!=NULL )
				color[i] = (float) rParam;
			else
				color[i] = 0;
		}
#endif
	}
	if( bProrieties[kStoppingVertex] && GeoVert.GetIsStoppingVertex() )
	{
		color[0] = color[2] = 0;
		color[1] = 1;
	}

#if 0	// just for debugging
	T_GeodesicVertexList& VertList = GW_Toolkit::GetStartVertex();
	if( !VertList.empty() )
	{
		GW_GeodesicVertex* pVert = VertList.front();
		GW_Float rDist = ~(pVert->GetPosition() - GeoVert.GetPosition());
		rDist = GW_ABS( rDist - GeoVert.GetDistance() );
		color[0] = color[1] = color[2] = rDist*5;
	}
#endif
}


/*------------------------------------------------------------------------------*/
// Name : GW_GeodesicDisplayer::DisplayPath
/**
 *  \param  CurPath [GW_GeodesicPath&] The path.
 *  \param  Color [GW_Vertex&] The color.
 *  \author Gabriel Peyré
 *  \date   4-10-2003
 * 
 *  Display a path.
 */
/*------------------------------------------------------------------------------*/
void GW_GeodesicDisplayer::DisplayPath( GW_GeodesicPath& CurPath, GW_Vector3D& Color, GW_Float rLineWidth )
{
	glLineWidth( (GLfloat) rLineWidth );
	glPointSize( (GLfloat) rLineWidth );
	T_GeodesicPointList& PointList = CurPath.GetPointList();
	glDisable( GL_LIGHTING );
	glColor( Color );
	glBegin( GL_LINE_STRIP );
	for( IT_GeodesicPointList it = PointList.begin(); it!=PointList.end(); ++it )
	{
		GW_GeodesicPoint* pPoint = *it;
		GW_ASSERT( pPoint->GetVertex1()!=NULL );
		GW_ASSERT( pPoint->GetVertex2()!=NULL );
		GW_Vector3D Pos = pPoint->GetVertex1()->GetPosition()*pPoint->GetCoord() + 
						  pPoint->GetVertex2()->GetPosition()*(1-pPoint->GetCoord());
		glColor( Color );
		glVertex( Pos );
		GW_GeodesicFace* pFace = pPoint->GetCurFace();
		GW_ASSERT( pFace!=NULL );
		T_SubPointVector& SubPointVector = pPoint->GetSubPointVector();
		GW_Vector3D& v0 = pPoint->GetVertex1()->GetPosition();
		GW_Vector3D& v1 = pPoint->GetVertex2()->GetPosition();
		GW_Vertex* pLastVert = pFace->GetVertex( *pPoint->GetVertex1(), *pPoint->GetVertex2() );
		GW_ASSERT( pLastVert!=NULL );
		GW_Vector3D& v2 = pLastVert->GetPosition();
		for( IT_SubPointVector it=SubPointVector.begin(); it!=SubPointVector.end(); ++it )
		{
			GW_Vector3D& coord = *it;
			Pos = v0*coord[0] + v1*coord[1] + v2*coord[2];
			glVertex( Pos );
		}
	}
	glEnd();
#if 0
	glColor3f( 0,0,0 );
	glPointSize(4);
	glBegin( GL_POINTS );
	for( IT_GeodesicPointList it = PointList.begin(); it!=PointList.end(); ++it )
	{
		GW_GeodesicPoint* pPoint = *it;
		GW_ASSERT( pPoint->GetVertex1()!=NULL );
		GW_ASSERT( pPoint->GetVertex2()!=NULL );
		GW_Vector3D Pos = pPoint->GetVertex1()->GetPosition()*pPoint->GetCoord() + 
			pPoint->GetVertex2()->GetPosition()*(1-pPoint->GetCoord());
		glVertex( Pos );
		GW_GeodesicFace* pFace = pPoint->GetCurFace();
		GW_ASSERT( pFace!=NULL );
		T_SubPointVector& SubPointVector = pPoint->GetSubPointVector();
		GW_Vector3D& v0 = pPoint->GetVertex1()->GetPosition();
		GW_Vector3D& v1 = pPoint->GetVertex2()->GetPosition();
		GW_Vertex* pLastVert = pFace->GetVertex( *pPoint->GetVertex1(), *pPoint->GetVertex2() );
		GW_ASSERT( pLastVert!=NULL );
		GW_Vector3D& v2 = pLastVert->GetPosition();
		for( IT_SubPointVector it=SubPointVector.begin(); it!=SubPointVector.end(); ++it )
		{
			GW_Vector3D& coord = *it;
			Pos = v0*coord[0] + v1*coord[1] + v2*coord[2];
			glVertex( Pos );
		}
	}
	glEnd();

	glEnable( GL_LIGHTING );
	glLineWidth( 1 );
#endif
}

/*------------------------------------------------------------------------------*/
// Name : GW_GeodesicDisplayer::DisplayGeodesicBoundaries
/**
 *  \param  VoronoiMesh [GW_VoronoiMesh&] The mesh.
 *  \param  Color [GW_Vector3D&] The color.
 *  \param  rLineWidth [GW_Float] Width of the lines.
 *  \author Gabriel Peyré
 *  \date   4-23-2003
 * 
 *  Display each boundaries of the voronoi mesh.
 */
/*------------------------------------------------------------------------------*/
void GW_GeodesicDisplayer::DisplayGeodesicBoundaries( GW_VoronoiMesh& VoronoiMesh, GW_Vector3D& Color, GW_Float rLineWidth )
{
	glLineWidth( (GLfloat) rLineWidth );
	T_VertexPathMap& BoundariesMap = VoronoiMesh.GetGeodesicBoundariesMap();
	for( IT_VertexPathMap it=BoundariesMap.begin(); it!=BoundariesMap.end(); ++it )
	{
		T_GeodesicVertexList* pList = it->second;
		GW_ASSERT( pList!=NULL );
		glDisable( GL_LIGHTING );
		glColor( Color );
		glBegin( GL_LINE_STRIP );
		for( IT_GeodesicVertexList VertIt = pList->begin(); VertIt!=pList->end(); ++VertIt )
		{
			GW_GeodesicVertex* pVert = *VertIt;
			GW_ASSERT( pVert!=NULL );
			glVertex( *pVert );
		}
		glEnd();
		glEnable( GL_LIGHTING );
	}
	glLineWidth( 1 );
}


/*------------------------------------------------------------------------------*/
// Name : GW_GeodesicDisplayer::DisplayVoronoiPoints
/**
 *  \param  VoronoiMesh [GW_VoronoiMeshBuilder&] The mesh builder.
 *  \author Gabriel Peyré
 *  \date   4-12-2003
 * 
 *  Display the set of points of the builder.
 */
/*------------------------------------------------------------------------------*/
void GW_GeodesicDisplayer::DisplayVoronoiPoints( GW_VoronoiMesh& VoronoiMesh, GW_Vector3D& Color, GW_Float rPointSize )
{
	GW_Bool bIsTextureOn = glIsEnabled( GL_TEXTURE_2D );
	glDisable( GL_TEXTURE_2D );

	T_GeodesicVertexList BaseVertexList = VoronoiMesh.GetBaseVertexList();
	glPointSize( (GLfloat) rPointSize );
	glColor( Color );
	glDisable( GL_LIGHTING );
	glBegin( GL_POINTS );
	for( IT_GeodesicVertexList it = BaseVertexList.begin(); it!=BaseVertexList.end(); ++it )
	{
		GW_GeodesicVertex* pVert = *it;
		glVertex( pVert->GetPosition() );
	}
	glEnd();
	glEnable( GL_LIGHTING );

	if( bIsTextureOn )
		glEnable( GL_TEXTURE_2D );
}



///////////////////////////////////////////////////////////////////////////////
//  Copyright (c) Gabriel Peyré
///////////////////////////////////////////////////////////////////////////////
//                               END OF FILE                                 //
///////////////////////////////////////////////////////////////////////////////

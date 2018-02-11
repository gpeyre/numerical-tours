
/*------------------------------------------------------------------------------*/
/** 
 *  \file   GW_BasicDisplayer.h
 *  \brief  Definition of class \c GW_BasicDisplayer
 *  \author Gabriel Peyré
 *  \date   3-28-2003
 */ 
/*------------------------------------------------------------------------------*/

#ifndef _GW_BASICDISPLAYER_H_
#define _GW_BASICDISPLAYER_H_

#include "../gw_core/GW_Config.h"
#include "../gw_core/GW_Mesh.h"
#include "GW_OpenGLHelper.h"
#include <GL/gl.h>


namespace GW {

/*------------------------------------------------------------------------------*/
/** 
 *  \class  GW_BasicDisplayer
 *  \brief  Basic display of a mesh. Can be overloaded.
 *  \author Gabriel Peyré
 *  \date   3-28-2003
 *
 *  Just flush faces to OpenGL.
 */ 
/*------------------------------------------------------------------------------*/

class GW_BasicDisplayer
{

public:

	GW_BasicDisplayer();
	virtual ~GW_BasicDisplayer();

	virtual void DisplayMesh(GW_Mesh& Mesh);
	virtual void DisplayFace(GW_Face& Face);
	virtual void ComputeColor( GW_Vertex& pVert, float* color );
	virtual void SetUpDraw(GW_Mesh& Mesh);

	#define GW_DISPLAYER_NBR_PTIES 17
	enum T_DisplayerProprieties
	{
		kNormal,
		kMinCurvDirection,
		kMaxCurvDirection,
		kMinCurv,
		kMaxCurv,
		kMeanCurv,
		kGaussianCurv,
		kMaxAbsCurv,
		kLighting,
		kGeodesicDistance,
		kBoundaries,
		kMarchingState,
		kForceMonoColor,
		kVertexParametrization,
		kStoppingVertex,
		kUseFlatLighting,
		kGeodesicDistanceStreamColor,
	};

	GW_Bool GetPropriety( T_DisplayerProprieties Pty );

	void SetDraw( T_DisplayerProprieties Pty, GW_Bool bVal );
	void EnableDraw( T_DisplayerProprieties Pty	);
	void DisableDraw( T_DisplayerProprieties Pty );
	void ToggleDraw( T_DisplayerProprieties Pty	);

	void SetVectorScaling(GW_Float rVectorScaling);
	GW_Float GetVectorScaling();

	void MultiplyContrast( T_DisplayerProprieties Pty, GW_Float Factor );
	void IncreaseContrast( T_DisplayerProprieties Pty, GW_Float Factor = 1.0/1.1f );
	void DecreaseContrast( T_DisplayerProprieties Pty, GW_Float Factor = 1.1f );

	typedef void (*T_ComputeColorCallback)( GW_Vertex& pVert, float* color );
	void RegisterComputeColorCallback(T_ComputeColorCallback pFunc);

    //-------------------------------------------------------------------------
    /** \name Vertex array management. */
    //-------------------------------------------------------------------------
    //@{
	void BuildVertexArray( GW_Mesh& Mesh );
	void BuildColorArray( GW_Mesh& Mesh );
    //@}


protected:

	/** the OpenGL interleaved array is GL_T2F_C4F_N3F_V3F */
	#define T_OFFSET 0	// texture coords offset
	#define C_OFFSET 2	// color offset
	#define N_OFFSET 6	// texture offset
	#define V_OFFSET 9	// vertex offset
	#define VCN_OFFSET 12	// total size

	GW_Bool bProrieties[GW_DISPLAYER_NBR_PTIES];

	GW_Float aMinValue_[GW_DISPLAYER_NBR_PTIES];
	GW_Float aMaxValue_[GW_DISPLAYER_NBR_PTIES];

	GW_Float rVectorScaling_;

	/** OpenGL vertex arrays management */
	GLfloat* VertexArray_;
	GLuint*  FaceArray_;
	GW_U32 nNbrVertex_;
	GW_U32 nNbrFace_;

	/** a callback for vertex color display */
	T_ComputeColorCallback pComputeColorCallback_;

};



} // End namespace GW


#ifdef GW_USE_INLINE
	#include "GW_BasicDisplayer.inl"
#endif

#endif // _GW_BASICDISPLAYER_H_


///////////////////////////////////////////////////////////////////////////////
//  Copyright (c) Gabriel Peyré
///////////////////////////////////////////////////////////////////////////////
//                               END OF FILE                                 //
///////////////////////////////////////////////////////////////////////////////

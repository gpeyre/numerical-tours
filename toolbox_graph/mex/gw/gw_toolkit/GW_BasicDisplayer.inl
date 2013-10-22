/*------------------------------------------------------------------------------*/
/** 
 *  \file   GW_BasicDisplayer.inl
 *  \brief  Inlined methods for \c GW_BasicDisplayer
 *  \author Gabriel Peyré
 *  \date   4-6-2003
 */ 
/*------------------------------------------------------------------------------*/

#include "GW_BasicDisplayer.h"

namespace GW {


/*------------------------------------------------------------------------------*/
// Name : GW_BasicDisplayer constructor
/**
*  \author Gabriel Peyré
*  \date   4-3-2003
* 
*  Constructor.
*/
/*------------------------------------------------------------------------------*/
GW_INLINE
GW_BasicDisplayer::GW_BasicDisplayer()
:	rVectorScaling_		( 0.1f ),
	VertexArray_		( NULL ),
	FaceArray_			( NULL ),
	pComputeColorCallback_	( NULL )
{
	for( GW_U32 i=0; i<GW_DISPLAYER_NBR_PTIES; ++i )
		bProrieties[i] = GW_False;
}

/*------------------------------------------------------------------------------*/
// Name : GW_BasicDisplayer destructor
/**
 *  \author Gabriel Peyré
 *  \date   4-27-2003
 * 
 *  Destructor.
 */
/*------------------------------------------------------------------------------*/
GW_INLINE
GW_BasicDisplayer::~GW_BasicDisplayer()
{
	GW_DELETEARRAY( VertexArray_ );
	GW_DELETEARRAY( FaceArray_ );
}



/*------------------------------------------------------------------------------*/
// Name : GW_BasicDisplayer::SetDraw
/**
*  \param  Pty [T_DisplayerProprieties] The propriety.
*  \param  bVal [GW_Bool] The value.
*  \author Gabriel Peyré
*  \date   4-3-2003
* 
*  Set wether to draw or not given propriety.
*/
/*------------------------------------------------------------------------------*/
GW_INLINE
void GW_BasicDisplayer::SetDraw( T_DisplayerProprieties Pty, GW_Bool bVal )
{
	bProrieties[ Pty ] = bVal;
}

/*------------------------------------------------------------------------------*/
// Name : GW_BasicDisplayer::EnableDraw
/**
*  \param  Pty [T_DisplayerProprieties] The propriety.
*  \param  bVal [GW_Bool] The value.
*  \author Gabriel Peyré
*  \date   4-3-2003
* 
*  Draw given propriety.
*/
/*------------------------------------------------------------------------------*/
GW_INLINE
void GW_BasicDisplayer::EnableDraw( T_DisplayerProprieties Pty )
{
	bProrieties[ Pty ] = GW_True;
}

/*------------------------------------------------------------------------------*/
// Name : GW_BasicDisplayer::DisableDraw
/**
*  \param  Pty [T_DisplayerProprieties] The propriety.
*  \param  bVal [GW_Bool] The value.
*  \author Gabriel Peyré
*  \date   4-3-2003
* 
*  Don't draw given propriety.
*/
/*------------------------------------------------------------------------------*/
GW_INLINE
void GW_BasicDisplayer::DisableDraw( T_DisplayerProprieties Pty  )
{
	bProrieties[ Pty ] = GW_False;
}

/*------------------------------------------------------------------------------*/
// Name : GW_BasicDisplayer::ToggleDraw
/**
*  \param  Pty [T_DisplayerProprieties] The propriety.
*  \param  bVal [GW_Bool] The value.
*  \author Gabriel Peyré
*  \date   4-3-2003
* 
*  Set wether to draw or not given propriety.
*/
/*------------------------------------------------------------------------------*/
GW_INLINE
void GW_BasicDisplayer::ToggleDraw( T_DisplayerProprieties Pty )
{
	bProrieties[ Pty ] = !bProrieties[ Pty ];
}


/*------------------------------------------------------------------------------*/
// Name : GW_BasicDisplayer::GetPropriety
/**
*  \param  Pty [T_DisplayerProprieties] The propriety
*  \return [GW_Bool] The state.
*  \author Gabriel Peyré
*  \date   4-3-2003
* 
*  Get the value of a propriety.
*/
/*------------------------------------------------------------------------------*/
GW_INLINE
GW_Bool GW_BasicDisplayer::GetPropriety( T_DisplayerProprieties Pty )
{
	return bProrieties[Pty];
}


/*------------------------------------------------------------------------------*/
// Name : GW_BasicDisplayer::SetVectorScaling
/**
*  \param  rVectorScaling [GW_Float] The new factor.
*  \author Gabriel Peyré
*  \date   4-6-2003
* 
*  Set the scaling factor for normal/curvature vectors.
*/
/*------------------------------------------------------------------------------*/
GW_INLINE
void GW_BasicDisplayer::SetVectorScaling(GW_Float rVectorScaling)
{
	rVectorScaling_ = rVectorScaling;
}

/*------------------------------------------------------------------------------*/
// Name : GW_BasicDisplayer::GetVectorScaling
/**
 *  \return [GW_Float] Factor.
 *  \author Gabriel Peyré
 *  \date   4-6-2003
 * 
 *  Get the factor used for display of normal/curvature vectors.
 */
/*------------------------------------------------------------------------------*/
GW_INLINE
GW_Float GW_BasicDisplayer::GetVectorScaling()
{
	return rVectorScaling_;
}


/*------------------------------------------------------------------------------*/
// Name : GW_BasicDisplayer::MultiplyContrast
/**
 *  \param  Pty [T_DisplayerProprieties] The propriety.
 *  \param  Factor [GW_Float] The scaling factor.
 *  \author Gabriel Peyré
 *  \date   4-6-2003
 * 
 *  Multiply the contrast of a given display propriety
 */
/*------------------------------------------------------------------------------*/
GW_INLINE
void GW_BasicDisplayer::MultiplyContrast( T_DisplayerProprieties Pty, GW_Float Factor )
{
	GW_Float m = (aMaxValue_[Pty] + aMinValue_[Pty])/2;
	GW_Float d = (aMaxValue_[Pty] - aMinValue_[Pty])/2;
	aMinValue_[Pty] = m - d*Factor;
	aMaxValue_[Pty] = m + d*Factor;
}


/*------------------------------------------------------------------------------*/
// Name : GW_BasicDisplayer::MultiplyContrast
/**
*  \param  Pty [T_DisplayerProprieties] The propriety.
*  \param  Factor [GW_Float] The scaling factor.
*  \author Gabriel Peyré
*  \date   4-6-2003
* 
*  Multiply the contrast of a given display propriety
*/
/*------------------------------------------------------------------------------*/
GW_INLINE
void GW_BasicDisplayer::IncreaseContrast( T_DisplayerProprieties Pty, GW_Float Factor )
{
	this->MultiplyContrast( Pty, Factor );
}


/*------------------------------------------------------------------------------*/
// Name : GW_BasicDisplayer::MultiplyContrast
/**
*  \param  Pty [T_DisplayerProprieties] The propriety.
*  \param  Factor [GW_Float] The scaling factor.
*  \author Gabriel Peyré
*  \date   4-6-2003
* 
*  Multiply the contrast of a given display propriety
*/
/*------------------------------------------------------------------------------*/
GW_INLINE
void GW_BasicDisplayer::DecreaseContrast( T_DisplayerProprieties Pty, GW_Float Factor )
{
	this->MultiplyContrast( Pty, Factor );
}

/*------------------------------------------------------------------------------*/
// Name : GW_BasicDisplayer::RegisterComputeColorCallback
/**
 *  \param  pFunc [T_ComputeColorCallback] User function.
 *  \author Gabriel Peyré
 *  \date   5-20-2003
 * 
 *  The user can override the color computations of the vertex.
 */
/*------------------------------------------------------------------------------*/
GW_INLINE
void GW_BasicDisplayer::RegisterComputeColorCallback(T_ComputeColorCallback pFunc)
{
	pComputeColorCallback_ = pFunc;
}


} // End namespace GW


///////////////////////////////////////////////////////////////////////////////
//  Copyright (c) Gabriel Peyré
///////////////////////////////////////////////////////////////////////////////
//                               END OF FILE                                 //
///////////////////////////////////////////////////////////////////////////////

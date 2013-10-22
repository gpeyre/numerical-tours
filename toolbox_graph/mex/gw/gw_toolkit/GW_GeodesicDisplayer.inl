/*------------------------------------------------------------------------------*/
/** 
 *  \file   GW_GeodesicDisplayer.inl
 *  \brief  Inlined methods for \c GW_GeodesicDisplayer
 *  \author Gabriel Peyré
 *  \date   4-10-2003
 */ 
/*------------------------------------------------------------------------------*/

#include "GW_GeodesicDisplayer.h"

namespace GW {

/*------------------------------------------------------------------------------*/
// Name : GW_GeodesicDisplayer::SetColorStreamFixedRadius
/**
 *  \param  val [GW_Float] <0 => scale the stream to the current distance map.
 *  \author Gabriel Peyré
 *  \date   9-29-2003
 * 
 *  Set the length of the stream.
 */
/*------------------------------------------------------------------------------*/
GW_INLINE
void GW_GeodesicDisplayer::SetColorStreamFixedRadius( GW_Float val )
{
	rColorStreamFixedRadius_ = val;
}



} // End namespace GW


///////////////////////////////////////////////////////////////////////////////
//  Copyright (c) Gabriel Peyré
///////////////////////////////////////////////////////////////////////////////
//                               END OF FILE                                 //
///////////////////////////////////////////////////////////////////////////////

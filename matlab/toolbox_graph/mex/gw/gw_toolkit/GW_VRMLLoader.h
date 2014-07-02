
/*------------------------------------------------------------------------------*/
/** 
 *  \file   GW_VRMLLoader.h
 *  \brief  Definition of class \c GW_VRMLLoader
 *  \author Gabriel Peyré
 *  \date   5-20-2003
 */ 
/*------------------------------------------------------------------------------*/

#ifndef _GW_VRMLLOADER_H_
#define _GW_VRMLLOADER_H_

#include "../gw_core/GW_Config.h"
#include "../gw_core/GW_Mesh.h"

namespace GW {

/*------------------------------------------------------------------------------*/
/** 
 *  \class  GW_VRMLLoader
 *  \brief  Load data from a VRML file.
 *  \author Gabriel Peyré
 *  \date   5-20-2003
 *
 *  Very simple loader.
 */ 
/*------------------------------------------------------------------------------*/

class GW_VRMLLoader
{

public:

	static GW_I32 Load(GW_Mesh& Mesh, const char *name, const char* mode = "rt", GW_U32 nExtraVertexPad = 2, GW_Bool bFlipFaces = GW_False );
	static GW_I32 Save( GW_Mesh& Mesh, const char *name );
private:

};

} // End namespace GW


#endif // _GW_VRMLLOADER_H_


///////////////////////////////////////////////////////////////////////////////
//  Copyright (c) Gabriel Peyré
///////////////////////////////////////////////////////////////////////////////
//                               END OF FILE                                 //
///////////////////////////////////////////////////////////////////////////////

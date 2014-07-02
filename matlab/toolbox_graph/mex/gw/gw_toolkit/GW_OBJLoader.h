/*------------------------------------------------------------------------------*/
/** 
 *  \file   GW_OBJLoader.h
 *  \brief  Definition of class \c GW_OBJLoader
 *  \author Gabriel Peyré
 *  \date   4-9-2003
 */ 
/*------------------------------------------------------------------------------*/

#ifndef __GW_OBJLoader__
#define __GW_OBJLoader__

#include "../gw_core/GW_Config.h"
#include "../gw_core/GW_Mesh.h"

GW_BEGIN_NAMESPACE

class GW_OBJLoader
{
public:

	static GW_I32 Load(GW_Mesh& Mesh, const char *name, const char* mode = "rt", GW_Bool bFlipFaces = GW_False );
	static GW_I32 Save( GW_Mesh& Mesh, const char *name );

};

GW_END_NAMESPACE

#endif	// #ifdef __GW_OBJLoader__



///////////////////////////////////////////////////////////////////////////////
//  Copyright (c) Gabriel Peyré
///////////////////////////////////////////////////////////////////////////////
//                               END OF FILE                                 //
///////////////////////////////////////////////////////////////////////////////

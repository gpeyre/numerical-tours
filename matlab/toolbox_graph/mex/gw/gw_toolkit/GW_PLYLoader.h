
/*------------------------------------------------------------------------------*/
/** 
 *  \file   GW_PLYLoader.h
 *  \brief  Definition of class \c GW_PLYLoader
 *  \author Gabriel Peyré
 *  \date   4-1-2003
 */ 
/*------------------------------------------------------------------------------*/

#ifndef _GW_PLYLOADER_H_
#define _GW_PLYLOADER_H_

#include "../gw_core/GW_Config.h"
#include "../gw_core/GW_Mesh.h"
extern "C" {
#include "ply/ply.h"
}
namespace GW {

/*------------------------------------------------------------------------------*/
/** 
 *  \class  GW_PLYLoader
 *  \brief  Loader for .ply polygon file.
 *  \author Gabriel Peyré
 *  \date   4-1-2003
 *
 *  PLY is a 3D standard from Stanford.
 */ 
/*------------------------------------------------------------------------------*/

class GW_PLYLoader
{

public:

	static GW_I32 Load(GW_Mesh& Mesh, const char *name, const char* mode = "rt", GW_U32 nExtraVertexPad = 2, GW_Bool bFlipFaces = GW_True );
	static GW_I32 Save( GW_Mesh& Mesh, const char *name, GW_Bool bAscii=GW_True );

private:

};

} // End namespace GW


#endif // _GW_PLYLOADER_H_


///////////////////////////////////////////////////////////////////////////////
//  Copyright (c) Gabriel Peyré
///////////////////////////////////////////////////////////////////////////////
//                               END OF FILE                                 //
///////////////////////////////////////////////////////////////////////////////

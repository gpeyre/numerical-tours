/*------------------------------------------------------------------------------*/
/** 
 *  \file   GW_OpenGLHelper.h
 *  \brief  Somme cool function for OpenGL linking
 *  \author Gabriel Peyré
 *  \date   5-30-2003
 */ 
/*------------------------------------------------------------------------------*/

#ifndef _GW_OPENGLHELPER_H_
#define _GW_OPENGLHELPER_H_

#include "../gw_core/GW_Config.h"
#include "../gw_core/GW_Vertex.h"
#include <GL/gl.h>

namespace GW {


/** some opengl cool linking features */
inline
void glVertex( GW::GW_VectorStatic<3,GW_Float> &v )
{
	glVertex3f( (GLfloat) v[0], (GLfloat) v[1], (GLfloat) v[2] );
}

inline
void glNormal( GW::GW_VectorStatic<3,GW_Float> &v )
{
	glNormal3f( (GLfloat) v[0], (GLfloat) v[1], (GLfloat) v[2] );
}


inline
void glColor( GW::GW_VectorStatic<3,GW_Float> &v )
{
	glColor3f( (GLfloat) v[0], (GLfloat) v[1], (GLfloat) v[2] );
}

inline
void glVertex( GW::GW_Vertex &vert )
{
	glVertex( vert.GetPosition() );
}

inline
void glNormal( GW::GW_Vertex &vert )
{
	glNormal( vert.GetNormal() );
}

inline
void glTexCoord( GW::GW_Vertex &vert )
{
	GLfloat u = (GLfloat) vert.GetTexCoordU();
	GLfloat v = (GLfloat) vert.GetTexCoordV();
	glTexCoord2f( u, v );
}



} // End namespace GW

#endif // _GW_OPENGLHELPER_H_


///////////////////////////////////////////////////////////////////////////////
//  Copyright (c) Gabriel Peyré
///////////////////////////////////////////////////////////////////////////////
//                               END OF FILE                                 //
///////////////////////////////////////////////////////////////////////////////

/*------------------------------------------------------------------------------*/
/** 
 *  \file   GW_BasicDisplayer.cpp
 *  \brief  Definition of class \c GW_BasicDisplayer
 *  \author Gabriel Peyré
 *  \date   3-28-2003
 */ 
/*------------------------------------------------------------------------------*/


#ifdef GW_SCCSID
    static const char* sccsid = "@(#) GW_BasicDisplayer.cpp(c) Gabriel Peyré2003";
#endif // GW_SCCSID

#include "stdafx.h"
#include "GW_BasicDisplayer.h"

#ifndef GW_USE_INLINE
	#include "GW_BasicDisplayer.inl"
#endif

using namespace GW;

/*------------------------------------------------------------------------------*/
// Name : GW_BasicDisplayer::DisplayMesh
/**
 *  \param  Mesh [GW_Mesh&] The mesh.
 *  \author Gabriel Peyré
 *  \date   3-28-2003
 * 
 *  Basic display of a whole mesh.
 */
/*------------------------------------------------------------------------------*/
#include <GL/glut.h>
void GW_BasicDisplayer::DisplayMesh(GW_Mesh& Mesh)
{
	GW_Face* pFace = NULL;
	GW_Vertex* pVert = NULL;
	GW_Vector3D* vect = NULL;
	GW_Vector3D TempVect; 
	if( bProrieties[kLighting] && !bProrieties[kForceMonoColor] )
		glEnable(GL_LIGHTING);
	else
		glDisable(GL_LIGHTING);
	glColorMaterial( GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE );
	glEnable( GL_COLOR_MATERIAL );
	if( !bProrieties[kUseFlatLighting] )
	{
		#define SPECULAR_COLOR 1
		GLfloat aSpecular[4] = {SPECULAR_COLOR,SPECULAR_COLOR,SPECULAR_COLOR,1};
		glMaterialfv(GL_FRONT_AND_BACK, GL_SPECULAR, (GLfloat*) aSpecular);
		glMaterialf(GL_FRONT_AND_BACK, GL_SHININESS, 100);
		if( VertexArray_!=NULL && FaceArray_!=NULL )
		{
			/* flush vertex array */
#if 1
			if( !bProrieties[kForceMonoColor] )
				glInterleavedArrays( GL_T2F_C4F_N3F_V3F, 0, VertexArray_ );
			else
			{
				glVertexPointer( 3, GL_FLOAT, VCN_OFFSET*sizeof(GLfloat), &VertexArray_[V_OFFSET] );
				glNormalPointer( GL_FLOAT, VCN_OFFSET*sizeof(GLfloat), &VertexArray_[N_OFFSET] );
				glTexCoordPointer( 2, GL_FLOAT, VCN_OFFSET*sizeof(GLfloat), &VertexArray_[T_OFFSET] );
				glEnableClientState( GL_VERTEX_ARRAY );
				glEnableClientState( GL_NORMAL_ARRAY );
				glEnableClientState( GL_TEXTURE_COORD_ARRAY );
				glDisableClientState( GL_COLOR_ARRAY );
			}
			glDrawElements( GL_TRIANGLES, 3*nNbrFace_, GL_UNSIGNED_INT, FaceArray_ );
#else // zap long edge

			glInterleavedArrays( GL_T2F_C4F_N3F_V3F, 0, VertexArray_ );
			for( GW_U32 i=0; i<Mesh.GetNbrFace(); ++i )
			{
				GW_Face* pFace = Mesh.GetFace(i);	GW_ASSERT( pFace!=NULL );
				GW_Vector3D n = pFace->ComputeNormal();
				GW_Float l = ~( pFace->GetVertex(0)->GetPosition()-pFace->GetVertex(1)->GetPosition() );
				l = GW_MAX( l, ~( pFace->GetVertex(0)->GetPosition()-pFace->GetVertex(2)->GetPosition() ) );
				l = GW_MAX( l, ~( pFace->GetVertex(1)->GetPosition()-pFace->GetVertex(2)->GetPosition() ) );
				if( l<1.2 )
					glDrawElements( GL_TRIANGLES, 3, GL_UNSIGNED_INT, FaceArray_+3*i );
			}

			// buggy code ???
	/*		glBegin(GL_TRIANGLES);
			for( GW_U32 i=0; i<Mesh.GetNbrFace(); ++i )
			{
				GW_Face* pFace = Mesh.GetFace(i);	GW_ASSERT( pFace!=NULL );
				GW_Vector3D n = pFace->ComputeNormal();
				GW_Float l = ~( pFace->GetVertex(0)->GetPosition()-pFace->GetVertex(1)->GetPosition() );
				l = GW_MAX( l, ~( pFace->GetVertex(0)->GetPosition()-pFace->GetVertex(2)->GetPosition() ) );
				l = GW_MAX( l, ~( pFace->GetVertex(1)->GetPosition()-pFace->GetVertex(2)->GetPosition() ) );
				if( l<1.2 )
				{
					float color[4];
					for( GW_U32 k=0; k<3; ++k )
					{
						ComputeColor( *pFace->GetVertex(k), (float*) color );
						glColor3fv( (float*) color);
						glNormal( *pFace->GetVertex(k) );
						glTexCoord2d( pFace->GetVertex(k)->GetTexCoordU(), pFace->GetVertex(k)->GetTexCoordV() );
						//					glTexCoord( *pFace->GetVertex(0) );
						glVertex( *pFace->GetVertex(k) );
					}
				}
			}
			glEnd(); */
#endif
		}
	}
	else
	{
		if( !bProrieties[kForceMonoColor] )
			glColor(GW_Vector3D(1,1,1));
		glBegin(GL_TRIANGLES);
		for( GW_U32 i=0; i<Mesh.GetNbrFace(); ++i )
		{
			GW_Face* pFace = Mesh.GetFace(i);	GW_ASSERT( pFace!=NULL );
			GW_Vector3D n = pFace->ComputeNormal();
			glNormal(n);
			glVertex( *pFace->GetVertex(0) );
			glNormal(n);
			glVertex( *pFace->GetVertex(1) );
			glNormal(n);
			glVertex( *pFace->GetVertex(2) );
		}
		glEnd();
	}

	/* display the edges */
	if( bProrieties[kBoundaries] )
	{
		GW_Bool bIsTextureOn = glIsEnabled( GL_TEXTURE_2D );
		glDisable( GL_TEXTURE_2D );

		glDisable( GL_LIGHTING );
		glColor3f( 1,0,0 );
		glLineWidth( 3 );
		glBegin( GL_LINES );
		for( GW_U32 i = 0; i<Mesh.GetNbrFace(); ++i )
		{
			pFace = Mesh.GetFace(i);
			GW_ASSERT(pFace!=NULL);
			for( GW_U32 s=0; s<3; ++s )
			{
				if( pFace->GetFaceNeighbor(s)==NULL )
				{	
					GW_Vertex* pVert1 = pFace->GetVertex((s+1)%3);	GW_ASSERT( pVert1!=NULL );
					GW_Vertex* pVert2 = pFace->GetVertex((s+2)%3);	GW_ASSERT( pVert2!=NULL );
					if( ~(pVert1->GetPosition()-pVert2->GetPosition())<0.2 )
					{
						glVertex( pVert1->GetPosition() );
						glVertex( pVert2->GetPosition() );
					}
				}
			}
		}
		glEnd();
		glEnable( GL_LIGHTING );
		if( bIsTextureOn )
			glEnable( GL_TEXTURE_2D );
	}

	/* draw the normals and so on */
	if( bProrieties[kNormal] || bProrieties[kMinCurvDirection] || bProrieties[kMaxCurvDirection] )
	{
		glDisable(GL_LIGHTING);
		glBegin( GL_LINES );
			for( GW_U32 i = 0; i<Mesh.GetNbrVertex(); ++i )
			{
				pVert = Mesh.GetVertex(i);
				GW_ASSERT(pVert!=NULL);
				if( pVert->GetFace()!=NULL )
				{
					if( bProrieties[kNormal] )
					{
						glColor3f( 0,1,0 );
						glVertex( pVert->GetPosition() );
						glVertex( pVert->GetPosition() + pVert->GetNormal()*rVectorScaling_ );
					}
					if( bProrieties[kMinCurvDirection] )
					{
						glColor3f( 0,0,1 );
						glVertex( pVert->GetPosition() );
						glVertex( pVert->GetPosition() + pVert->GetMinCurvDirection()*rVectorScaling_ );
					}
					if( bProrieties[kMaxCurvDirection] )
					{
						glColor3f( 1,0,0 );
						glVertex( pVert->GetPosition() );
						glVertex( pVert->GetPosition() + pVert->GetMaxCurvDirection()*rVectorScaling_ );
					}
				}
			}
		glEnd();
		glEnable(GL_LIGHTING);
	}
}

/*------------------------------------------------------------------------------*/
// Name : GW_BasicDisplayer::DisplayFace
/**
 *  \param  Face [GW_Face&] The face.
 *  \author Gabriel Peyré
 *  \date   3-28-2003
 * 
 *  Display a single face.
 */
/*------------------------------------------------------------------------------*/
void GW_BasicDisplayer::DisplayFace(GW_Face& Face)
{
	GW_Vertex* pVert = NULL;
	GW_Vector3D* vect = NULL;

	glBegin( GL_TRIANGLES );
	for( GW_U32 v=0; v<3; ++v )
	{
		pVert = Face.GetVertex(v);
		glNormal( *pVert );
		glVertex( *pVert );
	}
	glEnd();
}

/*------------------------------------------------------------------------------*/
// Name : GW_BasicDisplayer::ComputeColor
/**
 *  \param  pVert [GW_Vertex&] The vertex.
 *  \param  color [GLfloat*] The color.
 *  \author Gabriel Peyré
 *  \date   4-3-2003
 * 
 *  Assign a color to a vertex according to it's curv info.
 */
/*------------------------------------------------------------------------------*/
void GW_BasicDisplayer::ComputeColor( GW_Vertex& Vert, float* color )
{
	color[0] = color[1] = color[2] = 0;
	if( bProrieties[kMaxCurv] )
		color[0] = (float) GW_SCALE_01(Vert.GetMaxCurv(), aMinValue_[kMaxCurv], aMaxValue_[kMaxCurv]);
	if( bProrieties[kGaussianCurv] )
		color[1] = (float) GW_SCALE_01(Vert.GetGaussianCurv(), aMinValue_[kGaussianCurv], aMaxValue_[kGaussianCurv]);
	if( bProrieties[kMeanCurv] )
		color[1] = (float) GW_SCALE_01(Vert.GetMeanCurv(), aMinValue_[kMeanCurv], aMaxValue_[kMeanCurv]);
	if( bProrieties[kMinCurv] )
		color[2] = (float) 1-GW_SCALE_01(Vert.GetMinCurv(), aMinValue_[kMinCurv], aMaxValue_[kMinCurv]);
	if( bProrieties[kMaxAbsCurv] )
		color[0] = color[1] = color[2] = (float) GW_SCALE_01( GW_MAX(GW_ABS(Vert.GetMinCurv()),GW_ABS(Vert.GetMaxCurv())), aMinValue_[kMaxAbsCurv], aMaxValue_[kMaxAbsCurv]);
	if( !bProrieties[kMaxCurv]  && !bProrieties[kGaussianCurv] &&
	    !bProrieties[kMeanCurv] && !bProrieties[kMinCurv] && !bProrieties[kMaxAbsCurv] )
		color[0] = color[1] = color[2] = 1;
	if( pComputeColorCallback_ )
		pComputeColorCallback_(Vert, color);
}



/*------------------------------------------------------------------------------*/
// Name : GW_BasicDisplayer::SetUpDraw
/**
 *  \param  Mesh [GW_Mesh&] Mesh to draw.
 *  \author Gabriel Peyré
 *  \date   4-3-2003
 * 
 *  Compute the scaling factors for each propriety.
 */
/*------------------------------------------------------------------------------*/
void GW_BasicDisplayer::SetUpDraw(GW_Mesh& Mesh)
{
#if 1
	for( GW_U32 i=0; i<GW_DISPLAYER_NBR_PTIES; ++i )
	{
		aMinValue_[i] = 0;
		aMaxValue_[i] = 0;
	}
	/* compue mean in aMaxValue_ */
	for( GW_U32 i = 0; i<Mesh.GetNbrVertex(); ++i )
	{
		GW_Vertex* pVert = Mesh.GetVertex(i);
		aMaxValue_[kMinCurv] += pVert->GetMinCurv();
		aMaxValue_[kMaxCurv] += pVert->GetMaxCurv();
		aMaxValue_[kGaussianCurv] += pVert->GetGaussianCurv();
		aMaxValue_[kMeanCurv] += pVert->GetMeanCurv();
		aMaxValue_[kMaxAbsCurv] += pVert->GetMaxAbsCurv();
	}
	aMaxValue_[kMinCurv] /= Mesh.GetNbrVertex();
	aMaxValue_[kMaxCurv] /= Mesh.GetNbrVertex();
	aMaxValue_[kGaussianCurv] /= Mesh.GetNbrVertex();
	aMaxValue_[kMeanCurv] /= Mesh.GetNbrVertex();
	aMaxValue_[kMaxAbsCurv] /= Mesh.GetNbrVertex();

	/* comute mean deviation and store it in aMinValue_ */
	for( GW_U32 i = 0; i<Mesh.GetNbrVertex(); ++i )
	{
		GW_Vertex* pVert = Mesh.GetVertex(i);
		aMinValue_[kMinCurv] += (pVert->GetMinCurv()-aMaxValue_[kMinCurv])*(pVert->GetMinCurv()-aMaxValue_[kMinCurv]);
		aMinValue_[kMaxCurv] += (pVert->GetMaxCurv()-aMaxValue_[kMaxCurv])*(pVert->GetMaxCurv()-aMaxValue_[kMaxCurv]);
		aMinValue_[kGaussianCurv] += (pVert->GetGaussianCurv()-aMaxValue_[kGaussianCurv])*(pVert->GetGaussianCurv()-aMaxValue_[kGaussianCurv]);
		aMinValue_[kMeanCurv] += (pVert->GetMeanCurv()-aMaxValue_[kMeanCurv])*(pVert->GetMeanCurv()-aMaxValue_[kMeanCurv]);
		aMinValue_[kMaxAbsCurv] += (pVert->GetMaxAbsCurv()-aMaxValue_[kMaxAbsCurv])*(pVert->GetMaxAbsCurv()-aMaxValue_[kMaxAbsCurv]);
	}

	for( GW_U32 i=0; i<GW_DISPLAYER_NBR_PTIES; ++i )
	{
		GW_Float m = aMaxValue_[i];
		aMaxValue_[i] = m+sqrt( aMinValue_[i]/Mesh.GetNbrVertex() );
		aMinValue_[i] = m-sqrt( aMinValue_[i]/Mesh.GetNbrVertex() );
	}
#else
	/* comute mean deviation and store it in aMinValue_ */
	for( GW_U32 i=0; i<GW_DISPLAYER_NBR_PTIES; ++i )
	{
		aMinValue_[i] = 1e9;
		aMaxValue_[i] = -1e9;
	}
#define SETUP_PTY(name)		\
	if( pVert->Get##name()>aMinValue_[k##name] ) aMinValue_[k##name] = pVert->Get##name(); \
	if( pVert->Get##name()<aMaxValue_[k##name] ) aMaxValue_[k##name] = pVert->Get##name(); 

	for( GW_U32 i = 0; i<Mesh.GetNbrVertex(); ++i )
	{
		GW_Vertex* pVert = Mesh.GetVertex(i);
		SETUP_PTY(MinCurv);
		SETUP_PTY(MaxCurv);
		SETUP_PTY(GaussianCurv);
		SETUP_PTY(MeanCurv);
		SETUP_PTY(MaxAbsCurv);
	}
#endif
}



/*------------------------------------------------------------------------------*/
// Name : GW_BasicDisplayer::BuildColorArray
/**
 *  \author Gabriel Peyré
 *  \date   4-16-2003
 * 
 *  Compute the color for every vertex of the mesh.
 */
/*------------------------------------------------------------------------------*/
void GW_BasicDisplayer::BuildColorArray( GW_Mesh& Mesh )
{
	float color[3];
	GW_ASSERT( VertexArray_!=NULL );
	for( GW_U32 i=0; i<nNbrVertex_; ++i )
	{
		GW_Vertex* pVert = Mesh.GetVertex(i);
		GW_ASSERT( pVert!=NULL );
		ComputeColor( *pVert, color );
		VertexArray_[VCN_OFFSET*i+0+C_OFFSET] = color[0];
		VertexArray_[VCN_OFFSET*i+1+C_OFFSET] = color[1];
		VertexArray_[VCN_OFFSET*i+2+C_OFFSET] = color[2];
		VertexArray_[VCN_OFFSET*i+3+C_OFFSET] = 1;
	}
}

/*------------------------------------------------------------------------------*/
// Name : GW_BasicDisplayer::GW_BasicDisplayer::BuildVertexArray
/**
 *  \param  Mesh [GW_Mesh&] The mesh.
 *  \author Gabriel Peyré
 *  \date   4-16-2003
 * 
 *  Set up each vertex array of the mesh.
 */
/*------------------------------------------------------------------------------*/
void GW_BasicDisplayer::BuildVertexArray( GW_Mesh& Mesh )
{
	GW_DELETEARRAY( VertexArray_ );
	GW_DELETEARRAY( FaceArray_ );

	nNbrVertex_ = Mesh.GetNbrVertex();
	nNbrFace_ = Mesh.GetNbrFace();
	VertexArray_ = new GLfloat[nNbrVertex_*VCN_OFFSET];
	FaceArray_ = new GLuint[3*nNbrFace_];

	for( GW_U32 i=0; i<nNbrVertex_; ++i )
	{
		GW_Vertex* pVert = Mesh.GetVertex(i);
		GW_ASSERT( pVert!=NULL );
		VertexArray_[VCN_OFFSET*i+0+V_OFFSET] = (GLfloat) pVert->GetPosition()[0];
		VertexArray_[VCN_OFFSET*i+1+V_OFFSET] = (GLfloat) pVert->GetPosition()[1];
		VertexArray_[VCN_OFFSET*i+2+V_OFFSET] = (GLfloat) pVert->GetPosition()[2];
		VertexArray_[VCN_OFFSET*i+0+N_OFFSET] = (GLfloat) pVert->GetNormal()[0];
		VertexArray_[VCN_OFFSET*i+1+N_OFFSET] = (GLfloat) pVert->GetNormal()[1];
		VertexArray_[VCN_OFFSET*i+2+N_OFFSET] = (GLfloat) pVert->GetNormal()[2];
		VertexArray_[VCN_OFFSET*i+0+T_OFFSET] = (GLfloat) pVert->GetTexCoordU();
		VertexArray_[VCN_OFFSET*i+1+T_OFFSET] = (GLfloat) pVert->GetTexCoordV();
	}
	for( GW_U32 i=0; i<nNbrFace_; ++i )
	{
		GW_Face* pFace = Mesh.GetFace( i );
		GW_ASSERT( pFace!=NULL );
		for( GW_U32 j=0; j<3; ++j )
		{
			GW_Vertex* pVert = pFace->GetVertex(j);
			GW_ASSERT( pVert!=NULL );
			FaceArray_[3*i+j] = pVert->GetID();
		}
	}

	this->BuildColorArray( Mesh );
}



///////////////////////////////////////////////////////////////////////////////
//  Copyright (c) Gabriel Peyré
///////////////////////////////////////////////////////////////////////////////
//                               END OF FILE                                 //
///////////////////////////////////////////////////////////////////////////////

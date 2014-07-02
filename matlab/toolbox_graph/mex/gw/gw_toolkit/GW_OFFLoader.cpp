/*------------------------------------------------------------------------------*/
/** 
 *  \file   GW_OFFLoader.cpp
 *  \brief  Definition of class \c GW_OFFLoader
 *  \author Gabriel Peyré
 *  \date   5-20-2003
 */ 
/*------------------------------------------------------------------------------*/


#ifdef GW_SCCSID
    static const char* sccsid = "@(#) GW_OFFLoader.cpp(c) Gabriel Peyré2003";
#endif // GW_SCCSID

#include "stdafx.h"
#include "GW_OFFLoader.h"


using namespace GW;


/*------------------------------------------------------------------------------*/
// Name : GW_OFFLoader::Load
/**
*  \param  Mesh [GW_Mesh&] Mash to load data to.
*  \param  name [char*] File name.
*  \return [GW_I32] >0 : loading successful.
*  \author Gabriel Peyré
*  \date   4-1-2003
* 
*  Load data from a .PLY file.
*/
/*------------------------------------------------------------------------------*/
GW_I32 GW_OFFLoader::Load(GW_Mesh& Mesh, const char *name, const char* mode, GW_Bool bFlipFaces)
{
	FILE* pFile = NULL;
	pFile = fopen( name, mode );
	if( pFile ==NULL )
		return GW_Error_Opening_File;

	char data[255];
	GW_U32 nNbrVertex = 0;
	GW_U32 nNbrFace = 0;
	GW_U32 nTemp = 0;
	float Vect[3]; 
	GW_U32 nVertexNumber[3];
	GW_U32 nFaceValence = 0;
	GW_Vertex* pVert = NULL;
	GW_Face* pFace = NULL;

	T_FloatVector vX, vY, vZ;

	std::map<int, int> VertCorrespondance;
	GW_Bool bVertexDone = GW_False;

	while( !feof(pFile) )
	{
		fscanf(pFile, "%s", &data);

		if( strcmp(data, "#")==0 )
		{
			/* skip the whole line */
		}
		else if( strcmp(data, "OFF")==0 )
		{
			/* begin to read the data */
			fscanf(pFile, "%u %u %u", &nNbrVertex, &nNbrFace, &nTemp);
			Mesh.SetNbrVertex( nNbrVertex );
			Mesh.SetNbrFace( nNbrFace );
			for( GW_U32 i=0; i<nNbrVertex; ++i )
			{
				/* load the position */
				fscanf(pFile, "%f %f %f,", &Vect[0], &Vect[1], &Vect[2] );
				GW_Vector3D Pos(Vect[0], Vect[1], Vect[2]);
				pVert = &Mesh.CreateNewVertex();
				pVert->SetPosition( Pos );
				Mesh.SetVertex( i, pVert );
			}
			for( GW_U32 i=0; i<nNbrFace; ++i )
			{
				/* load the face */
				fscanf(pFile, "%u %u %u %u", &nTemp, &nVertexNumber[0], &nVertexNumber[1], &nVertexNumber[2] );
				GW_ASSERT( nTemp==3 );
				if( nTemp!=3 )
					return GW_ERROR;
				pFace = &Mesh.CreateNewFace();
				for( GW_U32 j=0; j<3; ++j )
				{
					pVert = Mesh.GetVertex( nVertexNumber[j] );
					GW_ASSERT( pVert!=NULL );
					if( !bFlipFaces )
						pFace->SetVertex( *pVert, j );
					else
						pFace->SetVertex( *pVert, 2-j );
				}
				Mesh.SetFace( i, pFace );
			}

			fclose(pFile);

			return GW_OK;
		}
	}

	return GW_ERROR;

}

/*------------------------------------------------------------------------------*/
// Name : GW_OFFLoader::Save
/**
 *  \param  Mesh [GW_Mesh&] The mesh.
 *  \param  name [char*] The file name.
 *  \return [GW_I32] Was the saving successful ?
 *  \author Gabriel Peyré
 *  \date   5-27-2003
 * 
 *  Save the mesh to a file.
 */
/*------------------------------------------------------------------------------*/
GW_I32 GW_OFFLoader::Save( GW_Mesh& Mesh, const char *name )
{
	FILE* pFile = NULL;
	pFile = fopen( name, "wt" );
	if( pFile ==NULL )
		return GW_Error_Opening_File;
#define WRITELN(str) fprintf( pFile, str "\n" )

	WRITELN("OFF");
	fprintf( pFile, "%u %u 0\n", Mesh.GetNbrVertex(), Mesh.GetNbrFace() );

	for( GW_U32 i=0; i<Mesh.GetNbrVertex(); ++i )
	{
		GW_Vertex* pVert = Mesh.GetVertex(i);
		GW_ASSERT( pVert!=NULL );
		GW_Vector3D& Pos = pVert->GetPosition();
		fprintf(pFile, "%f %f %f,\n", Pos[0], Pos[1], Pos[2] );

	}
	for( GW_U32 i=0; i<Mesh.GetNbrFace(); ++i )
	{
		GW_Face* pFace = Mesh.GetFace(i);
		GW_ASSERT( pFace!=NULL );
		fprintf(pFile, "3 %d %d %d\n", pFace->GetVertex(0)->GetID(), pFace->GetVertex(1)->GetID(), pFace->GetVertex(2)->GetID() );
	}
	fclose(pFile);

	return GW_OK;
}


///////////////////////////////////////////////////////////////////////////////
//  Copyright (c) Gabriel Peyré
///////////////////////////////////////////////////////////////////////////////
//                               END OF FILE                                 //
///////////////////////////////////////////////////////////////////////////////


/*------------------------------------------------------------------------------*/
/** 
 *  \file   GW_CSVLoader.h
 *  \brief  Definition of class \c GW_CSVLoader
 *  \author Gabriel Peyré
 *  \date   6-26-2003
 */ 
/*------------------------------------------------------------------------------*/

#ifndef _GW_CSVLOADER_H_
#define _GW_CSVLOADER_H_

#include "../gw_core/GW_Config.h"

namespace GW {

/*------------------------------------------------------------------------------*/
/** 
 *  \class  GW_CSVLoader
 *  \brief  A comma separating file loader.
 *  \author Gabriel Peyré
 *  \date   6-26-2003
 *
 *  Store the vertex and the faces in different files.
 */ 
/*------------------------------------------------------------------------------*/

class GW_CSVLoader
{

public:

	static GW_I32 LoadVertex( GW_Mesh& Mesh, const char *name )
	{
		FILE* pFile = NULL;

		pFile = fopen( "test.txt", "rt" );
		float Vect[3]; 
		fscanf(pFile, "%f", &Vect[0] );

		pFile = fopen( name, "rt" );
		if( pFile ==NULL )
			return GW_Error_Opening_File;

		GW_U32 nNbrVertex = 0;

		while( !feof(pFile) )
		{
			nNbrVertex++;
			fscanf(pFile, "%f,%f,%f", &Vect[0], &Vect[1], &Vect[2] );
		}
		fclose( pFile );
		pFile = fopen( name, "rt" );
		if( pFile ==NULL )
			return GW_Error_Opening_File;
		if( nNbrVertex!=Mesh.GetNbrVertex() )
			Mesh.SetNbrVertex( nNbrVertex );
		for( GW_U32 i=0; i<nNbrVertex; ++i )
		{
			GW_Vertex* pVert = Mesh.GetVertex(i);
			GW_ASSERT( pVert!=NULL );
			fscanf(pFile, "%f,%f,%f", &Vect[0], &Vect[1], &Vect[2] );
			pVert->SetPosition( GW_Vector3D((float) Vect[0], (float) Vect[1], (float) Vect[2]) );
		}

		fclose( pFile );

		return GW_OK;
	}

	static GW_I32 SaveVertex( GW_Mesh& Mesh, const char *name )
	{
		FILE* pFile = NULL;
		pFile = fopen( name, "wt" );
		if( pFile ==NULL )
			return GW_Error_Opening_File;

		for( GW_U32 i=0; i<Mesh.GetNbrVertex(); ++i )
		{
			GW_Vertex* pVert = Mesh.GetVertex(i);
			GW_ASSERT( pVert!=NULL );
			GW_Vector3D& Pos = pVert->GetPosition();
			fprintf(pFile, " %f,%f,%f", Pos[0], Pos[1], Pos[2] );
			if( i!=Mesh.GetNbrVertex() )
				fprintf(pFile, " \n");
		}

		fclose( pFile );

		return GW_OK;
	}

private:

};

} // End namespace GW


#endif // _GW_CSVLOADER_H_


///////////////////////////////////////////////////////////////////////////////
//  Copyright (c) Gabriel Peyré
///////////////////////////////////////////////////////////////////////////////
//                               END OF FILE                                 //
///////////////////////////////////////////////////////////////////////////////

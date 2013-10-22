/*------------------------------------------------------------------------------*/
/** 
 *  \file   GW_OBJLoader.cpp
 *  \brief  Definition of class \c GW_OBJLoader
 *  \author Gabriel Peyré
 *  \date   5-20-2003
 */ 
/*------------------------------------------------------------------------------*/


#ifdef GW_SCCSID
    static const char* sccsid = "@(#) GW_OBJLoader.cpp(c) Gabriel Peyré2003";
#endif // GW_SCCSID

#include "stdafx.h"
#include "GW_OBJLoader.h"


using namespace GW;


/*------------------------------------------------------------------------------*/
// Name : GW_OBJLoader::Load
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
GW_I32 GW_OBJLoader::Load(GW_Mesh& Mesh, const char *name, const char* mode, GW_Bool bFlipFaces)
{
	FILE* pFile = NULL;
	pFile = fopen( name, mode );
	if( pFile ==NULL )
		return GW_Error_Opening_File;

	char data[255];
	GW_U32 nTemp = 0;

	typedef std::list<GW_Vector3D> T_Vector3DList;
	typedef T_Vector3DList::iterator IT_Vector3DList;
	typedef GW_VectorStatic<3,GW_I32> T_U32Vector;
	typedef std::list<T_U32Vector> T_U32VectorList;
	typedef T_U32VectorList::iterator IT_U32VectorList;
	T_Vector3DList vertex_list;
	T_Vector3DList normal_list;
	T_Vector2DList texture_list;
	T_U32VectorList face_list;

	float Vect[3]; 

	GW_Bool bLoadingVertex = GW_False;
	GW_I32 nVertOffset = 0;
	GW_I32 nFutureOffset = 0;

	while( !feof(pFile) )
	{
		fscanf(pFile, "%s", &data);

		if( strcmp(data, "v")==0 )
		{
			if( !bLoadingVertex )
				nVertOffset = nFutureOffset;
			bLoadingVertex = GW_True;
			/* begin to read the data */
			fscanf(pFile, "%f %f %f", &Vect[0], &Vect[1], &Vect[2] );
			vertex_list.push_back( GW_Vector3D(Vect[0], Vect[1], Vect[2]) );
		}
		else
		{
			// record offset
			if( bLoadingVertex )
			{
				bLoadingVertex = GW_False;
				nFutureOffset = (GW_U32) vertex_list.size();
			}
		}


		if( strcmp(data, "vt")==0 )
		{
			/* begin to read the data */
			fscanf(pFile, "%f %f", &Vect[0], &Vect[1] );
			texture_list.push_back( GW_Vector2D(Vect[0], Vect[1]) );
		}
		else if( strcmp(data, "vn")==0 )
		{
			/* begin to read the data */
			fscanf(pFile, "%f %f %f", &Vect[0], &Vect[1], &Vect[2] );
			normal_list.push_back( GW_Vector3D(Vect[0], Vect[1], Vect[2]) );
		}
		else if( strcmp(data, "f")==0 )
		{
			T_U32Vector f;
			fscanf(pFile, "%u/%u/%u %u/%u/%u %u/%u/%u", &f[0], &nTemp, &nTemp, 
					&f[1], &nTemp, &nTemp, &f[2], &nTemp, &nTemp );
			/* f[0] += nVertOffset-1;
			f[1] += nVertOffset-1;
			f[2] += nVertOffset-1; */
			f[0] -= 1;
			f[1] -= 1;
			f[2] -= 1;
			if( f[0]>=0 )
				face_list.push_back( f );
		}
		else	// swap the line
		{
			fgets(data, 100, pFile);
		}
	}
	fclose(pFile);


	GW_U32 nNbrVertex = (GW_U32) vertex_list.size();
	GW_U32 nNbrFace = (GW_U32) face_list.size();
	GW_Bool bUseTexture = GW_False;
	if( texture_list.size()==nNbrVertex )
		bUseTexture = GW_True;
	GW_Bool bUseNormal = GW_False;
	if( normal_list.size()==nNbrVertex )
		bUseNormal = GW_True;

	Mesh.SetNbrVertex( nNbrVertex );
	Mesh.SetNbrFace( nNbrFace );

	// set up vertices
	IT_Vector3DList it_vertex = vertex_list.begin();
	IT_Vector2DList it_texture = texture_list.begin();
	IT_Vector3DList it_normal = normal_list.begin();
	GW_Vector2D texture;
	GW_Vector3D pos;
	GW_Vector3D normal;
	GW_U32 num = 0;
	while( it_vertex!=vertex_list.end() )
	{
		pos = *it_vertex;
		it_vertex++;
		if( bUseTexture )
		{
			texture = *it_texture;
			it_texture++;
		}
		if( bUseNormal )
		{
			normal = *it_normal;
			it_normal++;
		}
		GW_Vertex& vert = Mesh.CreateNewVertex();
		Mesh.SetVertex( num, &vert );
		vert.SetPosition(pos);
		vert.SetNormal(normal);
		vert.SetTexCoords( texture[0], texture[1] );
		num++;
	}
	num = 0;
	for( IT_U32VectorList it_face = face_list.begin(); it_face!=face_list.end(); ++it_face )
	{
		T_U32Vector f = *it_face;
		GW_Face& face = Mesh.CreateNewFace();
		Mesh.SetFace( num, &face );
		for( GW_U32 i=0; i<3; ++i )
		{
			GW_Vertex* pVert = Mesh.GetVertex(f[i]);	GW_ASSERT( pVert!=NULL );
			face.SetVertex( *pVert, i );
		}
		num++;
	}

	return GW_OK;
}

/*------------------------------------------------------------------------------*/
// Name : GW_OBJLoader::Save
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
GW_I32 GW_OBJLoader::Save( GW_Mesh& Mesh, const char *name )
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

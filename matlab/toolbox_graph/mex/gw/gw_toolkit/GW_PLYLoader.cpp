/*------------------------------------------------------------------------------*/
/** 
 *  \file   GW_PLYLoader.cpp
 *  \brief  Definition of class \c GW_PLYLoader
 *  \author Gabriel Peyré
 *  \date   4-1-2003
 */ 
/*------------------------------------------------------------------------------*/


#ifdef GW_SCCSID
    static const char* sccsid = "@(#) GW_PLYLoader.cpp(c) Gabriel Peyré 2003";
#endif // GW_SCCSID

#include "stdafx.h"
#include "GW_PLYLoader.h"

using namespace GW;

void big2littleendian(float& x)
{
	char* c = (char*) &x;
	char t[4] = { c[3],c[2],c[1],c[0] };
	x = (float) (*t);
}
void big2littleendian(int& x)
{
	char* c = (char*) &x;
	char t[4] = { c[3],c[2],c[1],c[0] };
	x = (int) (*t);
}

/*------------------------------------------------------------------------------*/
// Name : GW_PLYLoader::Load
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
GW_I32 GW_PLYLoader::Load(GW_Mesh& Mesh, const char *name, const char* mode, GW_U32 nExtraVertexPad, GW_Bool bFlipFaces)
{
	
	FILE* pFile = NULL;
	pFile = fopen( name, mode );
	if( pFile ==NULL )
		return GW_Error_Opening_File;

	/* pass the file to ply reader */
	int nelems, num_elems, nprops;
	char *elem_name;
	PlyFile* in_ply  = read_ply( pFile );
	if( in_ply==NULL ) 
		return GW_Error_Opening_File;
	char** elist = get_element_list_ply( in_ply, &nelems );

	/* are the data store using big endian ? */
	GW_Bool bConvertToLittleEndian = (in_ply->file_type==PLY_BINARY_BE);
	
	/* some temporary structures */
	struct plyVertex 
	{
		/* the usual 3-space position of a vertex */
		float x;
		float y;
		float z;
		/* the color */
		unsigned char red;
		unsigned char green;
		unsigned char blue;
	};

	struct plyFace 
	{
		unsigned char red;
		unsigned char green;
		unsigned char blue;
		unsigned char nverts;   // number of vertex indices in list
		int *verts;             // vertex index list
	};
	/* the proprieties we use to get back the data */
	PlyProperty vertProps[] = {
		{"x", Float32, Float32, static_cast<int>(offsetof(plyVertex,x)), 
			0, 0, 0, 0},
		{"y", Float32, Float32, static_cast<int>(offsetof(plyVertex,y)), 
			0, 0, 0, 0},
		{"z", Float32, Float32, static_cast<int>(offsetof(plyVertex,z)), 
			0, 0, 0, 0},
		{"red", Uint8, Uint8, static_cast<int>(offsetof(plyVertex,red)), 0, 0, 0, 0},
		{"green", Uint8, Uint8, static_cast<int>(offsetof(plyVertex,green)), 0, 0, 0, 0},
		{"blue", Uint8, Uint8, static_cast<int>(offsetof(plyVertex,blue)), 0, 0, 0, 0},
	};
	PlyProperty faceProps[] = {
		{"vertex_indices", Int32, Int32, 
			static_cast<int>(offsetof(plyFace,verts)),
			1, Uint8, Uint8, static_cast<int>(offsetof(plyFace,nverts))},
		{"red", Uint8, Uint8, static_cast<int>(offsetof(plyFace,red)), 0, 0, 0, 0},
		{"green", Uint8, Uint8, static_cast<int>(offsetof(plyFace,green)), 0, 0, 0, 0},
		{"blue", Uint8, Uint8, static_cast<int>(offsetof(plyFace,blue)), 0, 0, 0, 0},
	};

	/* Check to make sure that we can read geometry */
	PlyElement *elem;
	int index;
	if( (elem = find_element (in_ply, "vertex")) == NULL ||
		find_property (elem, "x", &index) == NULL ||
		find_property (elem, "y", &index) == NULL ||
		find_property (elem, "z", &index) == NULL ||
		(elem = find_element (in_ply, "face")) == NULL ||
			find_property (elem, "vertex_indices", &index) == NULL )
	{
		cerr << "Cannot read geometry." << endl;
		close_ply(in_ply);
		return GW_ERROR;
	}

	/*	Check for optional attribute data. We can handle the color red, green, blue. */
	GW_Bool bFaceColorAvailable = GW_False, bVertexColorAvailable = GW_False;
	if( (elem = find_element (in_ply, "face")) != NULL &&
		find_property (elem, "red", &index) != NULL &&
		find_property (elem, "green", &index) != NULL &&
		find_property (elem, "blue", &index) != NULL )
		bFaceColorAvailable = GW_True;

	if( (elem = find_element( in_ply, "vertex" )) != NULL &&
		find_property (elem, "red", &index) != NULL &&
		find_property (elem, "green", &index) != NULL &&
		find_property (elem, "blue", &index) != NULL )
		bVertexColorAvailable = GW_True;

	/* Okay, now we can grab the data */
	for( GW_I32 i = 0; i < nelems; ++i ) 
	{
		/* get the description of the first element */
		elem_name = elist[i];

		get_element_description_ply(in_ply, elem_name, &num_elems, &nprops);
// old		ply_get_element_description(in_ply, elem_name, &num_elems, &nprops);

		/* if we're on vertex elements, read them in */
		if( elem_name && strcmp("vertex", elem_name)==0 ) 
		{
			GW_U32 nNbrVertex = num_elems;
			/* set the number of vertex */
			Mesh.SetNbrVertex( nNbrVertex );
			
			ply_get_property(in_ply, elem_name, &vertProps[0]);
			ply_get_property(in_ply, elem_name, &vertProps[1]);
			ply_get_property(in_ply, elem_name, &vertProps[2]);

			/* Setup to read the PLY elements */
			if( bVertexColorAvailable )
			{
				ply_get_property(in_ply, elem_name, &vertProps[3]);
				ply_get_property(in_ply, elem_name, &vertProps[4]);
				ply_get_property(in_ply, elem_name, &vertProps[5]);
			}
			
			/* read the vertices */
			plyVertex vertex;
			for( GW_U32 j=0; j <nNbrVertex; ++j ) 
			{
				ply_get_element(in_ply, (void*) &vertex);
				GW_Vertex* pVert = &Mesh.CreateNewVertex();
				if( bConvertToLittleEndian )
				{
					big2littleendian(vertex.x);
					big2littleendian(vertex.y);
					big2littleendian(vertex.z);
				}
				pVert->SetPosition( GW_Vector3D(vertex.x, vertex.y, vertex.z) );
				Mesh.SetVertex( j, pVert );
				if ( bVertexColorAvailable )
				{
					// pVert->SetColor( GW_Vector3D(vertex.red, vertex.green, vertex.blue) );
				}
			}
		}	//if vertex
		else if( elem_name && !strcmp ("face", elem_name) ) 
		{
			GW_U32 nNbrFace = num_elems;
			Mesh.SetNbrFace( nNbrFace );
			
			/* Setup to read the PLY elements */
			ply_get_property(in_ply, elem_name, &faceProps[0]);
			if( bFaceColorAvailable )
			{
				ply_get_property(in_ply, elem_name, &faceProps[1]);
				ply_get_property(in_ply, elem_name, &faceProps[2]);
				ply_get_property(in_ply, elem_name, &faceProps[3]);
			}
			/* grab all the face elements */
			for( GW_U32 j=0; j < nNbrFace; j++) 
			{
				// grab and element from the file
				plyFace face;
				int verts[256];		// just 3 should be ok, but ...
				face.verts = verts;
				ply_get_element( in_ply, (void*) &face );
				// if( bConvertToLittleEndian )
				//	big2littleendian(face.nverts);
		/*		GW_ASSERT( face.nverts==3 ); 
				if( face.nverts!=3 ) 
				{
					cerr << "Face with more than 3 vertices are not supported." << endl;
					return GW_ERROR;
				}*/
				GW_Face* pFace = &Mesh.CreateNewFace();
				Mesh.SetFace( j, pFace );
				for( GW_U32 k=0; k<3; ++k)
				{
					if( bConvertToLittleEndian )
						big2littleendian(face.verts[k]);
					GW_ASSERT( face.verts[k]>=0 );
					GW_Vertex* pVert = Mesh.GetVertex( face.verts[k] );
					GW_ASSERT( pVert!=NULL ); 
					if( !bFlipFaces ) 
						pFace->SetVertex( *pVert, k );
					else
						pFace->SetVertex( *pVert, 2-k );
				}
				if( bFaceColorAvailable )
				{
					// pFace->SetColor( GW_Vector3D(face.red, face.green, face.blue) );
				}
			}
		}	//if face
		free(elist[i]); //allocated by ply_open_for_reading
    }	//for all elements of the PLY file
	free(elist); //allocated by ply_open_for_reading

	/* close the PLY files */
	close_ply( in_ply );
	fclose( pFile );

	return GW_OK;
}



/*------------------------------------------------------------------------------*/
// Name : GW_PLYLoader::Save
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
GW_I32 GW_PLYLoader::Save( GW_Mesh& Mesh, const char *name, GW_Bool bAscii )
{
	
	/* some temporary structures */
	struct plyVertex 
	{
		/* the usual 3-space position of a vertex */
		float x;
		float y;
		float z;
	};

	struct plyFace 
	{
		unsigned char nverts;   // number of vertex indices in list
		int *verts;             // vertex index list
	};
	/* the proprieties we use to get back the data */
	PlyProperty vertProps[] = {
		{"x", Float32, Float32, static_cast<int>(offsetof(plyVertex,x)), 
			0, 0, 0, 0},
		{"y", Float32, Float32, static_cast<int>(offsetof(plyVertex,y)), 
			0, 0, 0, 0},
		{"z", Float32, Float32, static_cast<int>(offsetof(plyVertex,z)), 
			0, 0, 0, 0},
	};
	PlyProperty faceProps[] = {
		{"vertex_indices", Int32, Int32, 
			static_cast<int>(offsetof(plyFace,verts)),
			1, Uint8, Uint8, static_cast<int>(offsetof(plyFace,nverts))},
	};

	FILE* pFile = NULL;
	pFile = fopen( name, "wt" );
	if( pFile ==NULL )
		return GW_Error_Opening_File;

	int file_type;
	if( bAscii )
		file_type = PLY_ASCII;
	else
		file_type = PLY_BINARY_LE;
	const char* elem_name[2] = {"vertex", "face"};
	PlyFile* out_ply  = write_ply( pFile, 2, (char**) elem_name, file_type );

	if( out_ply==NULL ) 
		return GW_Error_Opening_File;



	element_count_ply( out_ply, "vertex", Mesh.GetNbrVertex() );
	ply_describe_property( out_ply, "vertex", &vertProps[0]);
	ply_describe_property( out_ply, "vertex", &vertProps[1]);
	ply_describe_property( out_ply, "vertex", &vertProps[2]);

	
	element_count_ply( out_ply, "face", Mesh.GetNbrFace() );
	ply_describe_property( out_ply, "face", &faceProps[0]);

	/* write a comment and an object information field */
	append_comment_ply( out_ply, "Generated by GeoWave -- http://nikopol0.altern.org/geowave/" );
	append_obj_info_ply( out_ply, "GeoWave: v1.0" );

	/* complete the header */
	header_complete_ply( out_ply );

	/* set up and write the vertex elements */
	plyVertex vert;
	put_element_setup_ply( out_ply, "vertex" );

	for( GW_U32 i=0; i<Mesh.GetNbrVertex(); ++i )
	{
		GW_Vertex* pVert = Mesh.GetVertex(i);
		GW_ASSERT( pVert!=NULL );
		GW_Vector3D& Pos = pVert->GetPosition();
		vert.x = (float) Pos[0];
		vert.y = (float) Pos[1];
		vert.z = (float) Pos[2];
		put_element_ply( out_ply, (void*) &vert );
	}

	/* set up and write the face elements */
	plyFace face;
	int verts[3];
	face.verts = verts;
	put_element_setup_ply( out_ply, "face" );
	for( GW_U32 i=0; i<Mesh.GetNbrFace(); ++i )
	{
		GW_Face* pFace = Mesh.GetFace(i);
		GW_ASSERT( pFace!=NULL );
		for( GW_U32 j=0; j<3; ++j )
		{
			GW_ASSERT( pFace->GetVertex(j)!=NULL );
			face.verts[j] = (int) pFace->GetVertex(j)->GetID();
		}
		face.nverts = 3;
		put_element_ply( out_ply, (void*) &face );
	}

	put_element_ply( out_ply, (void*) &face );
	fclose(pFile);

	return GW_OK;
}




///////////////////////////////////////////////////////////////////////////////
//  Copyright (c) Gabriel Peyré
///////////////////////////////////////////////////////////////////////////////
//                               END OF FILE                                 //
///////////////////////////////////////////////////////////////////////////////

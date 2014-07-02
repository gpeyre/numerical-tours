/*------------------------------------------------------------------------------*/
/** 
 *  \file   GW_ASELoader.cpp
 *  \brief  Definition of class \c GW_ASELoader
 *  \author Gabriel Peyré
 *  \date   4-9-2003
 */ 
/*------------------------------------------------------------------------------*/


#ifdef GW_SCCSID
    static const char* sccsid = "@(#) GW_ASELoader.cpp(c) Gabriel Peyré 2003";
#endif // GW_SCCSID

#include "stdafx.h"
#include "GW_ASELoader.h"

using namespace GW;

const char *GW_ASELoader::mNumVertex  = "*MESH_NUMVERTEX";
const char *GW_ASELoader::mNumFaces   = "*MESH_NUMFACES";
const char *GW_ASELoader::mNumTVertex = "*MESH_NUMTVERTEX";
const char *GW_ASELoader::mNumTFaces  = "*MESH_NUMTVFACES";
const char *GW_ASELoader::mVertexList = "*MESH_VERTEX_LIST";
const char *GW_ASELoader::mVertex     = "*MESH_VERTEX";
const char *GW_ASELoader::mFaceList   = "*MESH_FACE_LIST";
const char *GW_ASELoader::mFace       = "*MESH_FACE";
const char *GW_ASELoader::mNormals    = "*MESH_NORMALS";
const char *GW_ASELoader::mFaceNormal = "*MESH_FACENORMAL";
const char *GW_ASELoader::mNVertex    = "*MESH_VERTEXNORMAL";
const char *GW_ASELoader::mTVert      = "*MESH_TVERT";
const char *GW_ASELoader::mTFace      = "*MESH_TFACE";
const char *GW_ASELoader::mTexture    = "*BITMAP";
const char *GW_ASELoader::mUTile      = "*UVW_U_TILING";
const char *GW_ASELoader::mVTile      = "*UVW_V_TILING";
const char *GW_ASELoader::mUOffset    = "*UVW_U_OFFSET";
const char *GW_ASELoader::mVOffset    = "*UVW_V_OFFSET";

void GW_ASELoader::Allocate()
{
	verts = new GW_Real32[numVerts * 3];
	faces = new GW_U32[numFaces * 3];

	if (want & kTexCoord)
	{
		texVerts = new GW_Real32[numTexVerts * 2];
		texFaces = new GW_U32[numTexFaces * 3];
	}
}

void GW_ASELoader::Delete()
{
	delete [] verts;
	delete [] texVerts;
	delete [] vertNorms;
	delete [] faces;
	delete [] texFaces;

	verts = texVerts = vertNorms = 0;
	faces = texFaces = 0; 
}

void GW_ASELoader::GetData(FILE *strm)
{
	char data[255];

	rewind(strm);

	while (!feof (strm))
	{
		fscanf(strm, "%s", &data);

		if (!strcmp(data, mVertex))
			GetVertex(strm);
		else if (want & kTexCoord && !strcmp(data, mTVert))
			GetTVertex (strm);
		else if (!strcmp(data, mFace))
			GetFace(strm);
		else if (want & kTexCoord && !strcmp(data, mTFace))
			GetTFace(strm);
		else if (!strcmp(data, mUTile))
			texUTile = GetFloatVal(strm);
		else if (!strcmp(data, mVTile))
			texVTile = GetFloatVal(strm);
		else if (!strcmp(data, mUOffset))
			texUOffset = GetFloatVal(strm);
		else if (!strcmp (data, mVOffset))
			texVOffset = GetFloatVal(strm);
		else fgets(data, 255, strm);
	}
}

void GW_ASELoader::GetInfo(FILE *strm)
{
	char data[255];

	rewind(strm);

	while (!feof(strm))
	{  
		fscanf(strm, "%s", &data);

		if (!strcmp(data, mNumVertex))
			fscanf(strm, "%d", &numVerts);
		else if (!strcmp(data, mNumFaces))
			fscanf(strm, "%d", &numFaces);
		else if (!strcmp(data, mNumTVertex))
			fscanf(strm, "%d", &numTexVerts);
		else if (!strcmp(data, mNumTFaces))
			fscanf(strm, "%d", &numTexFaces);
		else fgets(data, 255, strm);
	}	

	if (want & kTexCoord && !numTexVerts)
		want &= ~kTexCoord;
}

void GW_ASELoader::GetVertex(FILE *strm)
{
	GW_I32 index;

	fscanf(strm, "%d", &index);
	// swap y and z cause 3dsm likes too
	fscanf(strm, "%f %f %f", &verts[index * 3 + 0], &verts[index * 3 + 2],
		&verts[index * 3 + 1]);

	// in 3dsm negative z goes out of the screen, we want it to go in
	verts[3 * index + 2] = -verts[3 * index + 2];
}

void GW_ASELoader::GetFace(FILE *strm)
{
	GW_I32 index;

	fscanf(strm, "%d:", &index);

	fscanf(strm, "\tA:\t%d B:\t%d C:\t%d", &faces[3 * index + 0],
		&faces[3 * index + 1], &faces[3 * index + 2]); 
}

void GW_ASELoader::GetTFace(FILE *strm)
{
	GW_I32 index;

	fscanf(strm, "%d", &index);

	fscanf(strm, "%d %d %d", &texFaces[index * 3 + 0], &texFaces[3 * index + 1], 
		&texFaces[3 * index + 2]); 
}

void GW_ASELoader::GetTVertex(FILE *strm)
{
	GW_I32 index;

	// we ignore the z value, no idea what it does
	fscanf(strm, "%d", &index);
	fscanf(strm, "%f %f", &texVerts[2 * index + 0], &texVerts[2 * index +1]);

	if (texUTile)
		texVerts[2 * index + 0] *= texUTile;

	if (texVTile)
		texVerts[2 * index + 1] *= texVTile;
}

GW_Real32 GW_ASELoader::GetFloatVal(FILE *strm)
{
	GW_Real32 v;

	fscanf (strm, " %f", &v);

	return v;
}

void GW_ASELoader::Align()
{
	GW_U32 *pFace = faces;
	GW_U32 *pTex = texFaces;
	GW_Real32 *newVerts;
	GW_Real32 *newNorms;

	numVerts = numTexVerts;
	newVerts = new GW_Real32[numVerts * 3];

	if (want & kNormals)
		newNorms = new GW_Real32[numVerts * 3];

	for (GW_I32 cnt = 0; cnt < numFaces * 3; cnt++, pFace++, pTex++)
	{
		newVerts[3 * pTex[0] + 0] = verts[3 * pFace[0] + 0];
		newVerts[3 * pTex[0] + 1] = verts[3 * pFace[0] + 1];
		newVerts[3 * pTex[0] + 2] = verts[3 * pFace[0] + 2];

		if (want & kNormals)
		{
			newNorms[3 * pTex[0] + 0] = vertNorms[3 * pFace[0] + 0];
			newNorms[3 * pTex[0] + 1] = vertNorms[3 * pFace[0] + 1];
			newNorms[3 * pTex[0] + 2] = vertNorms[3 * pFace[0] + 2];
		}
	}

	if (faces)
		delete [] faces;
	if (verts)
		delete [] verts;

	faces = texFaces;
	verts = newVerts;
	texFaces = 0;

	if (want & kNormals && vertNorms)
	{
		GW_DELETEARRAY( vertNorms );
		vertNorms = newNorms;
	}
}

void GW_ASELoader::MakeNormals()
{
	GW_Vector3D *fNorms = NULL;
	GW_Real32 *vNorms = NULL;
	GW_Vector3D a, b, c, p, q;
	GW_Vector3D sum(0, 0, 0);
	GW_Real32 *pt1 = NULL, *pt2 = NULL, *pt3 = NULL;
	GW_U32 *pFace;
	GW_I32 cnt, pos;
	GW_I32 incident = 0;

	fNorms = new GW_Vector3D[numFaces];
	vNorms = new GW_Real32[3*numVerts];

	// find face normals for our model, used to find vert normals
	for (cnt = 0; cnt < numFaces; cnt++)
	{
		pt1 = verts + faces[cnt * 3 + 0] * 3;
		pt2 = verts + faces[cnt * 3 + 1] * 3;
		pt3 = verts + faces[cnt * 3 + 2] * 3;
		a.SetCoord(pt1[0], pt1[1], pt1[2]);
		b.SetCoord(pt2[0], pt2[1], pt2[2]);
		c.SetCoord(pt3[0], pt3[1], pt3[2]);

		p = b - a;
		q = c - a;

		fNorms[cnt] = p ^ q;
		fNorms[cnt].Normalize();
	}

	// find vertex normals
	for (cnt = 0; cnt < numVerts; cnt++)
	{
		for(pos = 0; pos < numFaces; pos++)
		{
			pFace = faces + pos * 3;

			if(pFace[0] == (GW_U32) cnt || pFace[1] == (GW_U32) cnt ||
				pFace[2] == (GW_U32) cnt)
			{
				incident++;
				sum += fNorms[pos];
			}
		}

		sum = - sum / ((GW_Real32) incident);
		sum.Normalize();
		vNorms[3*cnt]	= (GW_Real32) sum[0];
		vNorms[3*cnt+1] = (GW_Real32) sum[1];
		vNorms[3*cnt+2] = (GW_Real32) sum[2];
		sum.SetCoord(0,0,0);
		incident = 0;
	}

	GW_DELETEARRAY( fNorms );

	vertNorms = (GW_Real32 *) vNorms;
}

GW_ASELoader::GW_ASELoader()
{
	want = numVerts = numFaces = numTexFaces = numTexVerts = numIndex = 0;
	verts = texVerts = vertNorms = 0;
	faces = texFaces = 0;
	texUTile = texVTile = texUOffset = texVOffset = 0;
}

GW_ASELoader::~GW_ASELoader()
{
	Delete();
}

void GW_ASELoader::Release()
{
	Delete();
	want = numVerts = numFaces = numTexFaces = numTexVerts = numIndex = 0;
	texUTile = texVTile = texUOffset = texVOffset = 0;
}

GW_I32 GW_ASELoader::Load(GW_Mesh& Mesh, const char *name, GW_I32 bits)
{
	FILE *strm = fopen(name, "r");

	if (!strm)
		return GW_Error_Opening_File;

	Release();

	want = bits;

	GetInfo(strm);
	Allocate();
	GetData(strm);

	fclose(strm);

	if (want & kNormals)
		MakeNormals();

	if (want & kTexCoord && numTexVerts)
		Align();
	else if (want & kTexCoord && !numTexVerts)
		return 0;

	numIndex = numFaces * 3;

	/* retrieve information */
	Mesh.SetNbrVertex( this->GetNumVerts() );
	Mesh.SetNbrFace( this->GetNumFaces() );
	GW_U32* pFace = this->GetFaces();
	GW_Real32* pVert = this->GetVerts();
	GW_Real32* pNormal = this->GetVertNormals();
	GW_Real32* pTexture = this->GetTexture();
	GW_Vertex* pCurVert = NULL;
	GW_Vector3D Pos, Normal;

	std::map<int, int> VertCorrespondance;

	/* load vertex */
	for( GW_I32 i=0; i<this->GetNumVerts(); ++i )
	{
		Pos = GW_Vector3D( pVert[3*i], pVert[3*i+1], pVert[3*i+2] );
		/* this is to avoid an annoying problem of vertex duplication */
		for( GW_I32 j=0; j<i; ++j )
		{
			if( GW_ABS(pVert[3*i+0] - pVert[3*j+0])<GW_EPSILON && 
				GW_ABS(pVert[3*i+1] - pVert[3*j+1])<GW_EPSILON && 
				GW_ABS(pVert[3*i+2] - pVert[3*j+2])<GW_EPSILON )
				VertCorrespondance[j] = i;
		}
		if( pNormal!=NULL )
			Normal = GW_Vector3D( pNormal[3*i], pNormal[3*i+1], pNormal[3*i+2] );
		pCurVert = &Mesh.CreateNewVertex();
		pCurVert->SetPosition( Pos );
		pCurVert->SetNormal( Normal );
		if( want & kTexCoord )
		{
			for( GW_U32 s=0; s<1; ++s )
			{
				if( pTexture[2*i+s]<0 )
					pTexture[2*i+s] += 1;
				if( pTexture[2*i+s]>1 )
					pTexture[2*i+s] -= 1;
			}
			pCurVert->SetTexCoords( pTexture[2*i+0], pTexture[2*i+1] );
		} 
		Mesh.SetVertex( i, pCurVert );
	}
	/* load faces */
	GW_Face* pCurFace = NULL;
	for( GW_I32 i=0; i<this->GetNumFaces(); ++i )
	{
		pCurFace = &Mesh.CreateNewFace();
		GW_U32 FaceNumber[3];
		for( GW_U32 s=0; s<3; ++s  )
		{
			FaceNumber[s] = pFace[3*i+s];
			if( VertCorrespondance.find(FaceNumber[s])!=VertCorrespondance.end() )
				FaceNumber[s] = VertCorrespondance[FaceNumber[s]];
		}

		pCurFace->SetVertex( *Mesh.GetVertex(FaceNumber[0]), 
			*Mesh.GetVertex(FaceNumber[1]), 
			*Mesh.GetVertex(FaceNumber[2]) );
		Mesh.SetFace( i, pCurFace );
	}

	return GW_OK;
}



///////////////////////////////////////////////////////////////////////////////
//  Copyright (c) Gabriel Peyré
///////////////////////////////////////////////////////////////////////////////
//                               END OF FILE                                 //
///////////////////////////////////////////////////////////////////////////////

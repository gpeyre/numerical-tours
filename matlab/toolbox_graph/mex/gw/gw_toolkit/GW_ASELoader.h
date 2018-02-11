/*------------------------------------------------------------------------------*/
/** 
 *  \file   GW_ASELoader.h
 *  \brief  Definition of class \c GW_ASELoader
 *  \author Gabriel Peyré
 *  \date   4-9-2003
 */ 
/*------------------------------------------------------------------------------*/

#ifndef __GW_ASELoader__
#define __GW_ASELoader__

#include "../gw_core/GW_Config.h"
#include "../gw_core/GW_Mesh.h"

GW_BEGIN_NAMESPACE

class GW_ASELoader
{
public:

	GW_ASELoader();
	virtual ~GW_ASELoader();

	GW_I32 Load(GW_Mesh& Mesh, const char *name, GW_I32 bits = kTexCoord | kNormals);

protected:

	class aseFace_t
	{
	public:
		aseFace_t(GW_I32 n1, GW_I32 n2, GW_I32 n3) : ndx1(n1), ndx2(n2), ndx3(n3) {}
		GW_I32 ndx1, ndx2, ndx3;
	};


	void Release(); 

	GW_Real32 *GetVerts() const {return verts;}  
	GW_U32 *GetFaces() const {return faces;}
	GW_Real32 *GetVertNormals() const {return vertNorms;}
	GW_Real32 *GetTexture() const {return texVerts;}
	GW_I32 GetNumFaces() const {return numFaces;}
	GW_I32 GetNumVerts() const {return numVerts;}

	aseFace_t GetFace(GW_I32 num) const
	{
		return aseFace_t(3 * num + 0, 3 * num + 1, 3 * num + 2);
	}

	GW_Vector3D GetVertex(GW_I32 num) const
	{
		return GW_Vector3D(verts[3 * num + 0], verts[3 * num + 1], verts[3 * num + 2]);
	}

	GW_I32 GetNumIndex() {return numIndex;}

	enum T_LoadingFlag 
	{	
		kTexCoord = 0x01, 
		kNormals = 0x02
	}; 

	void Allocate();
	void Delete();
	void GetData (FILE *strm);
	void GetInfo(FILE *strm);
	void GetVertex(FILE *strm);
	void GetFace(FILE *strm);
	void GetTFace(FILE *strm);
	void GetTVertex (FILE *strm);
	GW_Real32 GetFloatVal (FILE *strm);
	void Align();
	void MakeNormals();

	static const char *mNumVertex;
	static const char *mNumFaces;
	static const char *mNumTVertex;
	static const char *mNumTFaces;
	static const char *mVertexList;
	static const char *mVertex;
	static const char *mFaceList;
	static const char *mFace;
	static const char *mNormals;
	static const char *mFaceNormal;
	static const char *mNVertex;
	static const char *mTVert;
	static const char *mTFace;
	static const char *mTexture;
	static const char *mUTile;
	static const char *mVTile;
	static const char *mUOffset;
	static const char *mVOffset;

	GW_I32 want;
	GW_I32 numVerts;
	GW_I32 numFaces;
	GW_I32 numTexFaces;
	GW_I32 numTexVerts;
	GW_I32 numIndex;
	GW_Real32 *verts;
	GW_Real32 *texVerts;
	GW_Real32 *vertNorms;
	GW_U32 *faces; 
	GW_U32 *texFaces; 

	GW_Real32 texUTile;
	GW_Real32 texVTile;
	GW_Real32 texUOffset;
	GW_Real32 texVOffset;
};

GW_END_NAMESPACE

#endif	// #ifdef __GW_ASELoader__



///////////////////////////////////////////////////////////////////////////////
//  Copyright (c) Gabriel Peyré
///////////////////////////////////////////////////////////////////////////////
//                               END OF FILE                                 //
///////////////////////////////////////////////////////////////////////////////

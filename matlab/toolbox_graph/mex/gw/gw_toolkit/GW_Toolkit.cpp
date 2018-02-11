/*------------------------------------------------------------------------------*/
/** 
 *  \file   GW_Toolkit.cpp
 *  \brief  Definition of class \c GW_Toolkit
 *  \author Gabriel Peyré
 *  \date   4-27-2003
 */ 
/*------------------------------------------------------------------------------*/


#ifdef GW_SCCSID
    static const char* sccsid = "@(#) GW_Toolkit.cpp(c) Gabriel Peyré2003";
#endif // GW_SCCSID

#include "stdafx.h"
#include "GW_Toolkit.h"

using namespace GW;

T_GeodesicVertexList GW_Toolkit::StartVertexList_;

const GW_Float GW_Toolkit::rRandomColor_[GW_TOOLKIT_NBR_COLOR][3] = 
{
	{1.0, 0.0, 0.0},
	{0.0, 1.0, 0.0},
	{0.0, 0.0, 1.0},
	{0.0, 1.0, 1.0},
	{1.0, 0.0, 1.0},
	{1.0, 1.0, 0.0},
	{1.0, 0.5, 0.5},
	{0.5, 1.0, 0.5},
	{0.5, 0.5, 1.0}, 


	{1.0, 0.5, 0.0}, 
	{1.0, 0.0, 0.5}, 
	{0.5, 1.0, 0.0}, 
	{0.0, 1.0, 0.5}, 
	{0.5, 0.0, 1.0}, 
	{0.0, 0.5, 1.0},

	{1.0, 1.0, 0.5},
	{1.0, 0.5, 1.0},
	{0.5, 1.0, 1.0}, 
};



/*------------------------------------------------------------------------------*/
// Name : GW_Toolkit constructor
/**
 *  \author Gabriel Peyré
 *  \date   4-27-2003
 * 
 *  Constructor.
 */
/*------------------------------------------------------------------------------*/
GW_Toolkit::GW_Toolkit()
:	nNbrStartVertex_( 1 )
{
	Displayer_.EnableDraw( GW_BasicDisplayer::kGeodesicDistance );
	Displayer_.SetVectorScaling( GW_TOOLKIT_VECTOR_SCALING );
	CoarseDisplayer_.EnableDraw( GW_BasicDisplayer::kLighting );
}

const GW_Bool bResavePLY = GW_False;
const GW_Bool bResaveVRML = GW_False;
const GW_Bool bResaveGWM = GW_True;
const GW_Bool bResaveOFF = GW_True;
const GW_Bool bResaveOBJ = GW_True;

/*------------------------------------------------------------------------------*/
// Name : GW_Toolkit::LoadFile
/**
*  \param  file_name [char*] File name.
*  \return [GW_I32] Was the loading successful ?
*  \author Gabriel Peyré
*  \date   4-2-2003
* 
*  Load a file depending on it's extension.
*/
/*------------------------------------------------------------------------------*/
GW_I32 GW_Toolkit::LoadFile(const char* file_name)
{
	StartVertexList_.clear();
	VoronoiMesh_.Reset();
	Mesh_.Reset();
	GeodesicPath_.ResetPath();
	pEndVertex_ = NULL;
	/* analyse extension */
	string str(file_name);
	string::size_type n = str.find_last_of( ".", str.length() );
	if( n!=string::npos )
	{
		string ext  = str.substr( n+1, str.size() );
		string base = str.substr( 0, n );
		if( ext=="ply" )
		{
#define MODE "rt"
#define EXTRA_PAD 8
#define FLIP_FACES GW_True
			GW_BasicDisplayer Displayer_;
			cout << "Loading PLY file " << file_name << "... ";
			GW_I32 nRet = GW_PLYLoader::Load( Mesh_, file_name, MODE, EXTRA_PAD, FLIP_FACES );
			if( nRet<0 )
			{
				cout << endl << "Can't load file.";
				return nRet;
			}
			cout << "done." << endl;
			/* set up the connectivity */
			cout << "Building connectivity ... ";
			Mesh_.BuildConnectivity();
			cout << "done." << endl;
			/* rescaling */
			cout << "Rescaling vertex ... ";
			Mesh_.TranslateVertex( -Mesh_.GetBarycenter() );
			Mesh_.ScaleVertex( 5/Mesh_.GetBoundingRadius() );
			cout << "done." << endl;
			/* building normal AFTER connectivity */
			cout << "Building curvature informations ... ";
			Mesh_.BuildCurvatureData();
			Mesh_.BuildRawNormal();	// todo : fix this.
			cout << "done." << endl;
			/* re-saving the file */
			if( bResaveGWM )
			{
				cout << "Saving file in .gwm format ... ";
				string base_gwm = base + ".gwm";
				GW_Ofstream file( base_gwm.c_str(), std::ofstream::out | std::ofstream::binary );
				file << Mesh_;
				file.close();
				cout << "done." << endl;
			}
			/* resave in vrml format */
			if( bResaveVRML )
			{
				cout << "Saving file in .wrl format ... ";
				string base_wrl = base + ".wrl";
				GW_VRMLLoader::Save( Mesh_, base_wrl.c_str() );
				cout << "done." << endl;
			}

		}
		else if( ext=="ase" )
		{	
			GW_ASELoader ASELoader;
			cout << "Loading ASE file " << file_name << " ... ";
			GW_I32 nRet = ASELoader.Load( Mesh_, file_name );
			if( nRet<0 )
			{
				cout << endl << "Can't load file.";
				return nRet;
			}
			cout << "done." << endl;
			/* set up the connectivity */
			cout << "Building connectivity ... ";
			Mesh_.BuildConnectivity();
			cout << "done." << endl;
			/* rescaling */
			cout << "Rescaling vertex ... ";
			Mesh_.TranslateVertex( -Mesh_.GetBarycenter() );
			Mesh_.ScaleVertex( 5/Mesh_.GetBoundingRadius() );
			cout << "done." << endl;
			/* building normal AFTER connectivity */
			cout << "Building curvature informations ... ";
			Mesh_.BuildCurvatureData();
			Mesh_.BuildRawNormal();	// todo : fix this.
			cout << "done." << endl;
			/* re-saving the file */
			if( bResaveGWM )
			{
				cout << "Saving file in .gwm format ... ";
				string base_gwm = base + ".gwm";
				GW_Ofstream file( base_gwm.c_str(), std::ofstream::out | std::ofstream::binary );
				file << Mesh_;
				file.close();
				cout << "done." << endl;
			}
			/* resave in OFF format */
			if( bResaveOFF )
			{
				cout << "Saving file in .off format ... ";
				string base_wrl = base + ".off";
				GW_OFFLoader::Save( Mesh_, base_wrl.c_str() );
				cout << "done." << endl;
			}
			/* resave in PLY format */
			if( bResavePLY )
			{
				cout << "Saving file in .ply format ... ";
				string base_wrl = base + ".ply";
				GW_PLYLoader::Save( Mesh_, base_wrl.c_str(), GW_False );
				cout << "done." << endl;
			}
		}
		else if( ext=="wrl" )
		{	
			cout << "Loading VRML file " << file_name << " ... ";
			GW_I32 nRet = GW_VRMLLoader::Load( Mesh_, file_name );
			if( nRet<0 )
			{
				cout << endl << "Can't load file.";
				return nRet;
			}
			cout << "done." << endl;
			/* set up the connectivity */
			cout << "Building connectivity ... ";
			Mesh_.BuildConnectivity();
			cout << "done." << endl;
			/* rescaling */
			cout << "Rescaling vertex ... ";
			Mesh_.TranslateVertex( -Mesh_.GetBarycenter() );
			GW_Float rRadius = Mesh_.GetBoundingRadius();
			GW_ASSERT( rRadius>0 );
			Mesh_.ScaleVertex( 5/(rRadius+0.001) );
			cout << "done." << endl;
			/* building normal AFTER connectivity */
			cout << "Building curvature informations ... ";
			Mesh_.BuildCurvatureData();
			Mesh_.BuildRawNormal();	// todo : fix this.
			cout << "done." << endl;
			/* re-saving the file */
			if( bResaveGWM )
			{
				cout << "Saving file in .gwm format ... ";
				string base_gwm = base + ".gwm";
				GW_Ofstream file( base_gwm.c_str(), std::ofstream::out | std::ofstream::binary );
				file << Mesh_;
				file.close();
				cout << "done." << endl;
			}
		}
		else if( ext=="off" )
		{	
			cout << "Loading OFF file " << file_name << " ... ";
			GW_I32 nRet = GW_OFFLoader::Load( Mesh_, file_name );
			if( nRet<0 )
			{
				cout << endl << "Can't load file.";
				return nRet;
			}
			cout << "done." << endl;
			/* set up the connectivity */
			cout << "Building connectivity ... ";
			GW_Vertex* pVert = Mesh_.GetVertex(0);
			Mesh_.GetVertex(10);
			Mesh_.BuildConnectivity();
			pVert = Mesh_.GetVertex(0);
			Mesh_.GetVertex(10);
			cout << "done." << endl;
			/* rescaling */
			cout << "Rescaling vertex ... ";
			Mesh_.TranslateVertex( -Mesh_.GetBarycenter() );
			pVert = Mesh_.GetVertex(0);
			Mesh_.GetVertex(10);
			GW_Float rRadius = Mesh_.GetBoundingRadius();
			GW_ASSERT( rRadius>0 );
			Mesh_.ScaleVertex( 5/(rRadius+0.001) );
			cout << "done." << endl;
			/* building normal AFTER connectivity */
			cout << "Building curvature informations ... ";
			Mesh_.BuildCurvatureData();
			Mesh_.BuildRawNormal();
			cout << "done." << endl;
			/* re-saving the file */
			if( bResaveOFF )
			{
				cout << "Saving file in .gwm format ... ";
				string base_gwm = base + ".gwm";
				GW_Ofstream file( base_gwm.c_str(), std::ofstream::out | std::ofstream::binary );
				file << Mesh_;
				file.close();
				cout << "done." << endl;
			}
		}

		else if( ext=="obj" )
		{	
			cout << "Loading OBJ file " << file_name << " ... ";
			GW_I32 nRet = GW_OBJLoader::Load( Mesh_, file_name );
			if( nRet<0 )
			{
				cout << endl << "Can't load file.";
				return nRet;
			}
			cout << "done." << endl;
			/* set up the connectivity */
			cout << "Building connectivity ... ";
			GW_Vertex* pVert = Mesh_.GetVertex(0);
			Mesh_.GetVertex(10);
			Mesh_.BuildConnectivity();
			pVert = Mesh_.GetVertex(0);
			Mesh_.GetVertex(10);
			cout << "done." << endl;
			/* rescaling */
			cout << "Rescaling vertex ... ";
			Mesh_.TranslateVertex( -Mesh_.GetBarycenter() );
			pVert = Mesh_.GetVertex(0);
			Mesh_.GetVertex(10);
			GW_Float rRadius = Mesh_.GetBoundingRadius();
			GW_ASSERT( rRadius>0 );
			Mesh_.ScaleVertex( 5/(rRadius+0.001) );
			cout << "done." << endl;
			/* building normal AFTER connectivity */
			cout << "Building curvature informations ... ";
			Mesh_.BuildCurvatureData();
			Mesh_.BuildRawNormal();
			cout << "done." << endl;
			/* re-saving the file */
			if( bResaveOBJ )
			{
				cout << "Saving file in .gwm format ... ";
				string base_gwm = base + ".gwm";
				GW_Ofstream file( base_gwm.c_str(), std::ofstream::out | std::ofstream::binary );
				file << Mesh_;
				file.close();
				cout << "done." << endl;
			}
		}
		else  if( ext=="gwm" )
		{
			cout << "Loading GWM file "  << file_name <<  " ... ";
			GW_Ifstream file( file_name, std::ios_base::in | std::ios_base::binary);
			file >> Mesh_;
			file.close();
			cout << "done." << endl;
			/* building normal */
			cout << "Building curvature informations ... ";
			Mesh_.BuildCurvatureData();
			Mesh_.BuildRawNormal();	// todo : fix this.
			cout << "done." << endl;
		}
		else
		{
			cerr << "load_file: " << file_name << " is not a supported file name." << endl;
			return GW_Error_File_Not_Supported;
		}
	}
	else
	{
		cerr << "load_file: " << file_name << " is not a supported file name." << endl;
		return GW_Error_File_Not_Supported;
	}

	char s[100];
	sprintf( s, "Number of vertices : %d", Mesh_.GetNbrVertex() );
	GW_OutputComment( s );
	sprintf( s, "Number of faces : %d", Mesh_.GetNbrFace() );
	GW_OutputComment( s );

	Displayer_.SetUpDraw( Mesh_ );
	Displayer_.BuildVertexArray( Mesh_ );
	this->SetUpMarching();

	return GW_OK;
}

void GW_Toolkit::SetUpMarching( GW_GeodesicVertex* pUserStart )
{
	StartVertexList_.clear();
	GW_GeodesicVertex* pStartVertex = NULL;
	Mesh_.ResetGeodesicMesh();
	/* find start vertex */
	for( GW_U32 i=0; i<nNbrStartVertex_; ++i )
	{
		pStartVertex = (GW_GeodesicVertex*) Mesh_.GetRandomVertex();
		GW_ASSERT( pStartVertex!=NULL );
		Mesh_.AddStartVertex( *pStartVertex );
		Displayer_.AddFrontColor( *pStartVertex, GW_Toolkit::GetRandomColor() );
		StartVertexList_.push_back( pStartVertex );
	}
	if( pUserStart!=NULL )
	{
		Mesh_.AddStartVertex( *pUserStart );
		Displayer_.AddFrontColor( *pUserStart, GW_Toolkit::GetRandomColor() );
		StartVertexList_.push_back( pUserStart );
	}

	Mesh_.SetUpFastMarching();
	Displayer_.EnableDraw( GW_BasicDisplayer::kMarchingState );
	Displayer_.EnableDraw( GW_BasicDisplayer::kLighting );
}

void GW_Toolkit::ComputePath( GW_Bool bRandomizeVertex )
{
	if( Mesh_.IsFastMarchingFinished() )
	{
		/* compute a path */
		if( bRandomizeVertex || pEndVertex_==NULL )
			pEndVertex_ = (GW_GeodesicVertex*) Mesh_.GetRandomVertex( GW_False );
		GW_ASSERT( pEndVertex_!=NULL );

		cout << "Computing geodesic path ... ";
		GeodesicPath_.ComputePath( *pEndVertex_, 5000 );
		cout << "done." << endl;


		Displayer_.BuildVertexArray( Mesh_ );
	}
	else
		cout << "You must finished Fast Marching first." << endl;
}

/*------------------------------------------------------------------------------*/
// Name : GW_Toolkit::AddFurthestPoint
/**
 *  \param  nNbr [GW_U32] Nbr points to add.
 *  \author Gabriel Peyré
 *  \date   4-27-2003
 * 
 *  Add furthest points to the coarse mesh.
 */
/*------------------------------------------------------------------------------*/
void GW_Toolkit::AddFurthestPoint( GW_U32 nNbr, GW_Bool bUseRandomStartVertex, GW_Bool bAssignColor, GW_U32 nPrintEach )
{
	/* add points to the list */
	GW_VoronoiMesh::AddFurthestPointsIterate( VoronoiMesh_.GetBaseVertexList(), Mesh_, nNbr, bUseRandomStartVertex );
	/* assign colors */
	Displayer_.ResetFrontColor();
	if( bAssignColor )
	for( IT_GeodesicVertexList it = VoronoiMesh_.GetBaseVertexList().begin(); it!=VoronoiMesh_.GetBaseVertexList().end(); ++it )
	{
		GW_GeodesicVertex* pVert = *it;
		GW_ASSERT( pVert!=NULL );
		Displayer_.AddFrontColor( *pVert, GW_Toolkit::GetRandomColor() );
	}	
	Displayer_.SetUpDraw( Mesh_ );
	Displayer_.DisableDraw( GW_BasicDisplayer::kMarchingState );
	Displayer_.BuildColorArray( Mesh_ );
}

/*------------------------------------------------------------------------------*/
// Name : GW_Toolkit::GetRandomColor
/**
 *  \return [GW_Vector3D] The color.
 *  \author Gabriel Peyré
 *  \date   4-29-2003
 * 
 *  Get a random color.
 */
/*------------------------------------------------------------------------------*/
GW_Vector3D GW_Toolkit::GetRandomColor()
{
	GW_U32 nNum = (GW_U32) floor( GW_RAND * (GW_TOOLKIT_NBR_COLOR+1) );
	GW_ASSERT( nNum<=GW_TOOLKIT_NBR_COLOR );
	return GW_Vector3D( rRandomColor_[nNum][0], rRandomColor_[nNum][1], rRandomColor_[nNum][2] );

}


/*------------------------------------------------------------------------------*/
// Name : GW_Toolkit::LoadFileList
/**
 *  \param  FileList [T_StringList&] The list to be feeded.
 *	\param  file [const char*] name of the file from which we read the list.
 *  \return [GW_I32] >=0 if success.
 *  \author Gabriel Peyré
 *  \date   5-5-2003
 * 
 *  Load a list of file.
 */
/*------------------------------------------------------------------------------*/
GW_I32 GW_Toolkit::LoadFileList( const char* file, T_StringList& FileList )
{
	/* load the file listing */
	FILE* pFile = fopen( file, "rt" );
	if( pFile==NULL )
		return -1;
	char data[255];
	while( !feof(pFile) )
	{
		fscanf(pFile, "%s", &data);
		if( data[0]=='#' )
		{
			// skip line
			char c = 0;
			while( !feof(pFile) && !(c=='\n') )
				c = fgetc( pFile );
		}
		else
			FileList.push_back( string(data) );
	}
	fclose( pFile );

	return 0;
}


/*------------------------------------------------------------------------------*/
// Name : GW_Toolkit::Display
/**
 *  \param  gc [GW_GeometryCell&] The cell
 *  \author Gabriel Peyré
 *  \date   2-4-2004
 * 
 *  Display a gometry cell.
 */
/*------------------------------------------------------------------------------*/
void GW_Toolkit::Display( GW_GeometryCell& gc, GW_Float boundary_thick, GW_Float half_boundary_thick, 
						 GW_Vector3D boundary_color, GW_Vector3D half_boundary_color, GW_Vector3D fill_color )
{
	glColor( fill_color );
	glLineWidth(1);
	for( GW_U32 j=0; j<(gc.GetWidth()-1); ++j )
	{
		glBegin(GL_QUAD_STRIP);
		for( GW_U32 i=0; i<gc.GetHeigth(); ++i )
		{
			glNormal( gc.GetNormal(i,j) );
			glVertex( gc.GetData(i,j) );
			glNormal( gc.GetNormal(i,j+1) );
			glVertex( gc.GetData(i,j+1) );
		}
		glEnd();
	}
	if( half_boundary_thick>0 )
	{	
		glColor( half_boundary_color );
		glLineWidth(half_boundary_thick);
		glBegin(GL_LINE_STRIP);
		for( GW_U32 i=0; i<gc.GetWidth(); ++i )
			glVertex( gc.GetData(i,0) );
		glEnd();
		glBegin(GL_LINE_STRIP);
		for( GW_U32 j=0; j<gc.GetHeigth(); ++j )
			glVertex( gc.GetData(0,j) );
		glEnd();
	}
	if( boundary_thick>0 )
	{
		glColor( boundary_color );
		glLineWidth(boundary_thick);
		glBegin(GL_LINE_STRIP);
		for( GW_I32 i=0; i<gc.GetWidth(); ++i )
			glVertex( gc.GetData(i,gc.GetHeigth()-1) );
		glEnd();
		glBegin(GL_LINE_STRIP);
		for( GW_U32 j=0; j<gc.GetHeigth(); ++j )
			glVertex( gc.GetData(gc.GetWidth()-1,j) );
		glEnd();
	}
	glColor( fill_color );
	glLineWidth(1);
}


/*------------------------------------------------------------------------------*/
// Name : GW_Toolkit::Display
/**
 *  \param  ga [GW_GeometryAtlas&] The atlas.
 *  \author Gabriel Peyré
 *  \date   2-4-2004
 * 
 *  Display a whole atlas of cells.
 */
/*------------------------------------------------------------------------------*/
void GW_Toolkit::Display( GW_GeometryAtlas& ga, GW_Float boundary_thick, GW_Float half_boundary_thick, 
						 GW_Vector3D boundary_color, GW_Vector3D half_boundary_color, GW_Vector3D fill_color )
{
	T_GeometryCellVector& CellVector = ga.GetCellVector();
	for( IT_GeometryCellVector it = CellVector.begin(); it!=CellVector.end(); ++it )
	{
		GW_GeometryCell* pCell = *it;	GW_ASSERT( pCell!=NULL );
		GW_Toolkit::Display(*pCell, boundary_thick, half_boundary_thick, boundary_color, half_boundary_color, fill_color);
	}
}


///////////////////////////////////////////////////////////////////////////////
//  Copyright (c) Gabriel Peyré
///////////////////////////////////////////////////////////////////////////////
//                               END OF FILE                                 //
///////////////////////////////////////////////////////////////////////////////

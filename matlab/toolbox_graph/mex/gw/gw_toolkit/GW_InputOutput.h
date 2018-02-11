/*------------------------------------------------------------------------------*/
/** 
 *  \file   GW_InputOutput.h
 *  \brief  Definition of class \c GW_InputOutput
 *  \author Gabriel Peyré
 *  \date   6-25-2003
 */ 
/*------------------------------------------------------------------------------*/

#ifndef _GW_INPUTOUTPUT_H_
#define _GW_INPUTOUTPUT_H_

#include "../gw_core/GW_Config.h"
#include "trackball.h"
#include <GL/glut.h>

#ifdef GW_DEBUG
	#pragma comment(lib, "glut32.lib")
#else
	#pragma comment(lib, "glut32.lib")
#endif // GW_DEBUG


namespace GW {

/*------------------------------------------------------------------------------*/
/** 
 *  \class  GW_InputOutput
 *  \brief  A minimalistic helper that use \c GLUT.
 *  \author Gabriel Peyré
 *  \date   6-25-2003
 *
 *  Just handle mouse.
 */ 
/*------------------------------------------------------------------------------*/

class GW_InputOutput
{

public:

	typedef void (*T_IdleFunc)();
	static void DefaultIdle()
	{
		if( bSpinning )
			add_quats(lastquat, curquat, curquat);
		glutPostRedisplay();
	}

	static void SetIdleCallback( T_IdleFunc NewFunc )
	{
		IdleCallback_ = NewFunc;
	}

	static void Reshape( int w, int h )
	{
		if( h!=0 )
		{
			glViewport( 0, 0, w, h  );
			glMatrixMode( GL_PROJECTION );
			glLoadIdentity();
			gluPerspective(30,w/h,1.5,50000);
			glMatrixMode( GL_MODELVIEW );
			glLoadIdentity();
			rWidth = w;
			rHeight = h;
			rAspectRatio = (GLdouble) rWidth/(GLdouble) rHeight;
		}
	}

	/*------------------------------------------------------------------------------*/
	// Name : mouse
	/**
	*  \param  x [GW::GW_I32] x position of the pointer.
	*  \param  y [GW::GW_I32] y position of the pointer.
	*  \author Gabriel Peyré
	*  \date   11-24-2002
	* 
	*  To handle mouse movements.
	*/
	/*------------------------------------------------------------------------------*/
	static void MouseMotion(int x, int y)
	{
		if( bTranslating )
		{
			GW_Float rTranslationX = -rZoom*TRANS_FACTOR*( x-nMouseX );
			GW_Float rTranslationY = rZoom*TRANS_FACTOR*( y-nMouseY );
			GW_InputOutput::MoveTarget( rTranslationX, rTranslationY );
		}
		else if( bMoving ) 
		{
			trackball(lastquat,
				ROT_SPEED*(2.0 * nBeginX - rWidth) / rWidth,
				ROT_SPEED*(rHeight - 2.0 * nBeginY) / rHeight,
				ROT_SPEED*(2.0 * x - rWidth) / rWidth,
				ROT_SPEED*(rHeight - 2.0 * y) / rHeight
				);
			nBeginX = x;
			nBeginY = y;
			bSpinning = 1;
			glutIdleFunc(IdleCallback_);
		}
		else if( bZooming )
		{
			rZoom += ZOOM_FACTOR*( y-nMouseY );
			if( rZoom<=0 )
				rZoom=0;
			if( rZoom>=ZOOM_MAX )
				rZoom=ZOOM_MAX;
		}
	}

	static void MouseClick( int button , int state , int x , int y )
	{
		bTranslating = GW_False;
		if( button == GLUT_RIGHT_BUTTON && state == GLUT_DOWN )
			bRightButton_ = GW_True;
		if( button == GLUT_RIGHT_BUTTON && state == GLUT_UP ) 
			bRightButton_ = GW_False;
		if( button == GLUT_LEFT_BUTTON && state == GLUT_DOWN )
			bLeftButton_ = GW_True;
		if( button == GLUT_LEFT_BUTTON && state == GLUT_UP )
			bLeftButton_ = GW_False;

		if( bLeftButton_ && bRightButton_ )
			bTranslating = GW_True;
		else if( bRightButton_ ) 
		{
			bTranslating = GW_False;
			bSpinning = 0;
			trackball(lastquat, 0, 0, 0, 0);
			glutIdleFunc(NULL);
			bMoving = 1;
			nBeginX = x;
			nBeginY = y;
		} 
		else if( !bRightButton_ && !bLeftButton_) 
		{
			bTranslating = GW_False;
			bMoving = 0;
		} 
		else if( bLeftButton_ ) 
		{
			bTranslating = GW_False;
			bZooming = GW_True;
		}
		else if( !bLeftButton_ ) 
		{
			bTranslating = GW_False;
			bZooming = GW_False;
		}
		nMouseY = y;
		nMouseX = x;
	}

	static void MousePassive( int x , int y )
	{
		nMouseY = y;
		nMouseX = x;
	}


	static void Init( int argc, char *argv[], GW_U32 nWindowsSize = 100, const char* windowsname = "GW -- Test application" )
	{
		srand( (unsigned) time( NULL ) );
		/* init the fast sqrt table */
		makeInverseSqrtLookupTable();
		/* init the track ball */
		trackball(curquat, -0.5, 0.0, 0.0, 0.0);

		glutInit( &argc , argv );
		glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB | GLUT_DEPTH);
		glutInitWindowSize( nWindowsSize, nWindowsSize ) ;
		glutInitWindowPosition( 100, 100 ) ;
		glutCreateWindow( windowsname );

		glutMotionFunc( MouseMotion );
		glutMouseFunc( MouseClick );
		glutIdleFunc( IdleCallback_ );
		glutReshapeFunc( Reshape );
	
		/* init opengl */
		glEnable( GL_DEPTH_TEST );
		glEnable( GL_CULL_FACE );
		glCullFace( GL_BACK );

		/* load light */
		GLfloat light_color[4] = {0.3,0.3,0.3,1};
		glLightfv( GL_LIGHT0, GL_AMBIENT, light_color );
		glLightfv( GL_LIGHT0, GL_DIFFUSE, light_color );
		glLightfv( GL_LIGHT0, GL_SPECULAR, light_color );
		glEnable( GL_LIGHTING );
		glEnable( GL_LIGHT0 );
		glDisable( GL_NORMALIZE );

		/* for line display */
		glEnable(GL_POLYGON_OFFSET_FILL);
		glPolygonOffset(3,3);

		glClearColor( BACKGROUND_COLOR,BACKGROUND_COLOR,BACKGROUND_COLOR,BACKGROUND_COLOR );

	}

	static void PrepareDisplay()
	{

		glClear( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT );
		glLoadIdentity();

		GLfloat m[4][4];
		GLfloat invm[4][4];
		//	TB_Vector objectSpaceLightPosition;
		TB_Vector objectSpaceEyePosition;

		/* Given the current rotation quaternion, generate the matching rotation
		matrix and combine it with the modelview matrix. */
		build_rotmatrix(m, curquat);

		/* Because we know that "m" is a rotation matrix, we could just
		perform a tranpose here, but for an arbitrary matrix, we 
		need a full matrix inversion. */
		invertMatrixf(&invm[0][0], &m[0][0]);

		/* Transform the light position in world space into object space
		using the inverse of the object space to world space transform. */
		transformPosition(&objectSpaceEyePosition, &eyePos, invm);

		/* Load a persepctive projection normally. */
		glMatrixMode(GL_PROJECTION);
		glLoadIdentity();
		gluPerspective( /* field of view in degree */ 40.0,
			/* aspect ratio */ rAspectRatio,
			/* Z near */ 1.0, /* Z far */ 500.0f);

		GLdouble eyeLookAt[16];
		/* Compute a "look at" 4x4 matrix with the "eye" at the eye position. */
		buildLookAtMatrix(eyePos.x, eyePos.y, eyePos.z,  0,0,0,  0,1,0,  eyeLookAt);
		glMatrixMode(GL_MODELVIEW);
		glTranslatef( -TargetPosition_[0], -TargetPosition_[1], -TargetPosition_[2] );
		glMultMatrixd(eyeLookAt);

		/* posit light */
		glTranslatef(0,0,-rZoom);
		glMultMatrixf(&m[0][0]);
		GLfloat light_position[4] = { 8,8,0,1 };
		glLightfv( GL_LIGHT0 , GL_POSITION, light_position ) ;
	}

	static void TranslateTarget( GW_Vector3D& tr )
	{
		TargetPosition_ += tr;
	}
	static void MoveTarget( GW_Float x, GW_Float y )	// moving in the "look at" plane
	{
		TargetPosition_[0] += x;
		TargetPosition_[1] += y;
	}

#define CAM_FILE "param.svg"
	static void SaveCameraSettings()
	{
		FILE* pFile = fopen( CAM_FILE, "wt" );
		if( pFile==NULL )
			return;
		fprintf( pFile, "cur_quat=%f,%f,%f,%f\n", curquat[0], curquat[1], curquat[2], curquat[3] );
		fprintf( pFile, "eye_pos=%f,%f,%f\n", eyePos.x, eyePos.y, eyePos.z );
		fprintf( pFile, "zoom=%f\n", (float) rZoom );
		fclose( pFile );
	}
	static void LoadCameraSettings()
	{
		FILE* pFile = fopen( CAM_FILE, "rt" );
		if( pFile==NULL )
			return;
		fscanf( pFile, "cur_quat=%f,%f,%f,%f\n", & curquat[0], & curquat[1], & curquat[2], & curquat[3] );
		fscanf( pFile, "eye_pos=%f,%f,%f\n", & eyePos.x, & eyePos.y, & eyePos.z );
		float zoom;
		fscanf( pFile, "zoom=%f\n", &zoom );
		rZoom = zoom;
		fclose( pFile );
	}

	static void RotateView(GW_Float phi, float a[3])
	{
		float q[4];
		axis_to_quat(a, phi, q);
		add_quats(q, curquat, curquat);
	}

private:

	/** for the interface */
	static GW_Bool bLeftButton_, bRightButton_;
	static GW_Bool bSpinning, bMoving;
	static GW_I32 nBeginX, nBeginY;
	static GW_I32 rWidth, rHeight;
	static GLdouble rAspectRatio;
	static float curquat[4];
	static float lastquat[4];
	static TB_Vector eyePos;
	static GW_Float rZoom;
	static GW_Bool bZooming, bTranslating;
	static GW_I32 nMouseY, nMouseX;
	static GW_U32 ZOOM_MAX;
	static GW_Float TRANS_FACTOR;
	static GW_Float ZOOM_FACTOR;
	static GW_Float ROT_SPEED;
	static GW_Float BACKGROUND_COLOR;
	static GW_Vector3D TargetPosition_;

	static T_IdleFunc IdleCallback_;

};

} // End namespace GW


#endif // _GW_INOUTOUTPUT_H_


///////////////////////////////////////////////////////////////////////////////
//  Copyright (c) Gabriel Peyré
///////////////////////////////////////////////////////////////////////////////
//                               END OF FILE                                 //
///////////////////////////////////////////////////////////////////////////////

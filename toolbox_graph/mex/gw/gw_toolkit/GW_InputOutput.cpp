/*------------------------------------------------------------------------------*/
/** 
 *  \file   GW_InputOutput.cpp
 *  \brief  Definition of class \c GW_InputOutput
 *  \author Gabriel Peyré
 *  \date   6-25-2003
 */ 
/*------------------------------------------------------------------------------*/


#ifdef GW_SCCSID
    const char* sccsid = "@(#) GW_InputOutput.cpp(c) Gabriel Peyré2003";
#endif // GW_SCCSID

#include "stdafx.h"
#include "GW_InputOutput.h"

using namespace GW;


GW_Bool GW_InputOutput::bSpinning = GW_False;
GW_Bool GW_InputOutput::bMoving = GW_False;
GW_I32 GW_InputOutput::nBeginX;
GW_I32 GW_InputOutput::nBeginY;
GW_I32 GW_InputOutput::rWidth = 300;
GW_I32 GW_InputOutput::rHeight = 300;
GLdouble GW_InputOutput::rAspectRatio;
float GW_InputOutput::curquat[4];
float GW_InputOutput::lastquat[4];
TB_Vector GW_InputOutput::eyePos = { 0.0, 0.0, 3.0 };
GW_Float GW_InputOutput::rZoom = 10;
GW_Bool GW_InputOutput::bZooming = GW_False;
GW_Bool GW_InputOutput::bTranslating = GW_False;
GW_I32 GW_InputOutput::nMouseY = 0;
GW_I32 GW_InputOutput::nMouseX = 0;
GW_U32 GW_InputOutput::ZOOM_MAX = 2000;
GW_Float GW_InputOutput::ZOOM_FACTOR = 0.2f;
GW_Float GW_InputOutput::ROT_SPEED = 0.3f;
GW_Float GW_InputOutput::BACKGROUND_COLOR = 1;
GW_Vector3D GW_InputOutput::TargetPosition_;
GW_Bool GW_InputOutput::bLeftButton_ = GW_False;
GW_Bool GW_InputOutput::bRightButton_ = GW_False;
GW_Float GW_InputOutput::TRANS_FACTOR = 0.00006;
GW_InputOutput::T_IdleFunc GW_InputOutput::IdleCallback_ = GW_InputOutput::DefaultIdle;

///////////////////////////////////////////////////////////////////////////////
//  Copyright (c) Gabriel Peyré
///////////////////////////////////////////////////////////////////////////////
//                               END OF FILE                                 //
///////////////////////////////////////////////////////////////////////////////

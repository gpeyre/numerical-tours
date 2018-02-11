/*------------------------------------------------------------------------------*/
/** 
 *  \file   GW_Config.cpp
 *  \brief  Definition of class \c GW_Config
 *  \author Gabriel Peyré
 *  \date   5-1-2003
 */ 
/*------------------------------------------------------------------------------*/


#ifdef GW_SCCSID
    static const char* sccsid = "@(#) GW_Config.cpp(c) Gabriel Peyré2003";
#endif // GW_SCCSID

#include "stdafx.h"
#include "GW_Config.h"

using namespace GW;

string GW::gw_alinea = "  * ";
string GW::gw_endl = "\n";
FILE* GW::GW_OutputStream = NULL;

void GW::GW_OutputComment( const char* str )
{
	if( GW_OutputStream!=NULL )
	{
		string str2 = gw_alinea + str + gw_endl;
		fprintf( GW_OutputStream, str2.c_str() );
	}
}

///////////////////////////////////////////////////////////////////////////////
//  Copyright (c) Gabriel Peyré
///////////////////////////////////////////////////////////////////////////////
//                               END OF FILE                                 //
///////////////////////////////////////////////////////////////////////////////

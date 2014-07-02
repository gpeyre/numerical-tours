/*------------------------------------------------------------------------------*/
/** 
 *  \file   GW_FaceIterator.cpp
 *  \brief  Definition of class \c GW_FaceIterator
 *  \author Gabriel Peyré
 *  \date   4-1-2003
 */ 
/*------------------------------------------------------------------------------*/


#ifdef GW_SCCSID
    static const char* sccsid = "@(#) GW_FaceIterator.cpp(c) Gabriel Peyré2003";
#endif // GW_SCCSID

#include "stdafx.h"
#include "GW_FaceIterator.h"
#include "GW_Face.h"

using namespace GW;

GW_FaceIterator::GW_FaceIterator(  GW_Face* pFace, GW_Vertex* pOrigin, GW_Vertex* pDirection, GW_U32 nNbrIncrement )
:	pFace_		( pFace), 
	pOrigin_	( pOrigin  ), 
	pDirection_	( pDirection ),
	nNbrIncrement_	( nNbrIncrement )
{ }

/* assignement */
GW_FaceIterator& GW_FaceIterator::operator=( const GW_FaceIterator& it)
{
	this->pFace_		= it.pFace_;
	this->pOrigin_		= it.pOrigin_;
	this->pDirection_	= it.pDirection_;
	this->nNbrIncrement_= it.nNbrIncrement_;
	return *this;
}

/* egality */
GW_Bool GW_FaceIterator::operator==( const GW_FaceIterator& it)
{
	return (this->pFace_==it.pFace_)
		&& (this->pOrigin_==it.pOrigin_)
		&& (this->pDirection_==it.pDirection_);
}

/* egality */
GW_Bool GW_FaceIterator::operator!=( const GW_FaceIterator& it)
{
	return (this->pFace_!=it.pFace_)
		|| (this->pOrigin_!=it.pOrigin_)
		|| (this->pDirection_!=it.pDirection_);
}


/* dereference */
GW_Face* GW_FaceIterator::operator*(  )
{
	return this->pFace_;
}

/* progression : \todo take in acount NULL pointer */
void GW_FaceIterator::operator++()
{
	if( nNbrIncrement_>100 )
	{
		GW_ASSERT( GW_False );
		(*this) = GW_FaceIterator(NULL,NULL,NULL);
		return;
	}

	if( pFace_!=NULL && pDirection_!=NULL && pOrigin_!=NULL )
	{
		GW_Face* pNextFace = pFace_->GetFaceNeighbor( *pDirection_ );
		/* check for end() */
		if(  pNextFace==pOrigin_->GetFace() )
		{
			(*this) = GW_FaceIterator(NULL,NULL,NULL);
		}
		else
		{
			if( pNextFace==NULL )
			{
				/* we are on a border face : Rewind on the first face */
				GW_Face* pPrevFace = pFace_;
				pDirection_ = pFace_->GetVertex( *pDirection_, *pOrigin_ );	// get rewind direction
				GW_ASSERT( pDirection_!=NULL );

				GW_U32 nIter = 0;
				do 
				{					
					pFace_ = pPrevFace;
					pPrevFace = pPrevFace->GetFaceNeighbor( *pDirection_ );
					pDirection_ = pFace_->GetVertex( *pOrigin_, *pDirection_ ); // next direction
					nIter++;
					GW_ASSERT( nIter<20 );
					if( nIter>=20 )
					{
						// this is on non-manifold ...
						(*this) = GW_FaceIterator(NULL,NULL,NULL);
						return;
					}

				}
				while( pPrevFace!=NULL );

				if( pFace_==pOrigin_->GetFace() )
				{
					// we are on End.
					(*this) = GW_FaceIterator(NULL,NULL,NULL);
				}
				else
				{
					GW_ASSERT( pDirection_!=NULL );
					(*this) = GW_FaceIterator( pFace_, pOrigin_, pDirection_, nNbrIncrement_+1 );
				}
				return;
			}
			GW_Vertex* pNextDirection = pFace_->GetVertex( *pOrigin_, *pDirection_ );
			GW_ASSERT( pNextDirection!=NULL );
			(*this) = GW_FaceIterator( pNextFace, pOrigin_, pNextDirection, nNbrIncrement_+1 );
		}
	}
	else
	{
		(*this) = GW_FaceIterator(NULL,NULL,NULL);
	}
}

/*------------------------------------------------------------------------------*/
// Name : GW_FaceIterator::GetLeftVertex
/**
 *  \return [GW_Vertex*] The vertex.
 *  \author Gabriel Peyré
 *  \date   5-6-2003
 * 
 *  Get the vertex on the left.
 */
/*------------------------------------------------------------------------------*/
GW_Vertex* GW_FaceIterator::GetLeftVertex()
{
	return pDirection_;
}

/*------------------------------------------------------------------------------*/
// Name : GW_FaceIterator::GetRightVertex
/**
*  \return [GW_Vertex*] The vertex.
*  \author Gabriel Peyré
*  \date   5-6-2003
* 
*  Get the vertex on the right.
*/
/*------------------------------------------------------------------------------*/
GW_Vertex* GW_FaceIterator::GetRightVertex()
{
	if( pFace_==NULL )
		return NULL;
	return pFace_->GetVertex( *pDirection_, *pOrigin_ );
}


///////////////////////////////////////////////////////////////////////////////
//  Copyright (c) Gabriel Peyré
///////////////////////////////////////////////////////////////////////////////
//                               END OF FILE                                 //
///////////////////////////////////////////////////////////////////////////////

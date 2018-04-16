/*****************************************************************************
* This file is provided under the Creative Commons Attribution 3.0 license.
*
* You are free to share, copy, distribute, transmit, or adapt this work
* PROVIDED THAT you attribute the work to the authors listed below.
* For more information, please see the following web page:
* http://creativecommons.org/licenses/by/3.0/
*
* This file is a component of the Sleipnir library for functional genomics,
* authored by:
* Curtis Huttenhower (chuttenh@princeton.edu)
* Mark Schroeder
* Maria D. Chikina
* Olga G. Troyanskaya (ogt@princeton.edu, primary contact)
*
* If you use this library, the included executable tools, or any related
* code in your work, please cite the following publication:
* Curtis Huttenhower, Mark Schroeder, Maria D. Chikina, and
* Olga G. Troyanskaya.
* "The Sleipnir library for computational functional genomics"
*****************************************************************************/
#ifndef HALFMATRIXI_H
#define HALFMATRIXI_H

#undef int64_t
#include <stdint.h>

namespace Sleipnir {

class CHalfMatrixBase {
protected:
	CHalfMatrixBase( ) : m_iSize(0) { }

	static void HalfIndex( size_t& iX, size_t& iY ) {
		size_t	i;

		if( iX > iY ) {
			i = iX;
			iX = iY;
			iY = i - iY - 1; }
		else
			iY -= iX + 1; }

	size_t	m_iSize;
};

}

#endif // HALFMATRIXI_H

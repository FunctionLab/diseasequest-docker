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
#ifndef STDAFX_H
#define STDAFX_H

#define __STDC_LIMIT_MACROS

#include <fstream>
#include <cstdlib>
#include <iomanip>
#include <queue>
#include <sstream>
using namespace std;

#ifdef _MSC_VER
#include <io.h>
#include <winsock2.h>
#else // _MSC_VER
#include <arpa/inet.h>
#include <netinet/in.h>

#include <omp.h>

#define SOCKET		int

inline bool _mktemp_s( char* szTemplate ) {

	return !mktemp( szTemplate ); }
#endif // _MSC_VER

#include <pthread.h>

#include <boost/graph/graphviz.hpp>
#undef INTMAX_C
#undef UINTMAX_C

#include "annotation.h"
#include "bayesnet.h"
#include "color.h"
#include "database.h"
#include "genome.h"
#include "meta.h"
#include "server.h"
#include "serverclient.h"
#include "statistics.h"
using namespace Sleipnir;

#endif // STDAFX_H

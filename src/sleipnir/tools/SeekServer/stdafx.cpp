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
#include "stdafx.h"

/*!
 * \page SeekServer SeekServer
 * 
 * SeekServer runs the coexpression mining algorithm using a multithreaded TCP/IP interface.
 * When it is running, SeekServer services requests from multiple clients over the network. The requests that can be
 * handled by SeekServer are, for example:
 * \li retrieve a list of genes that are found to be coexpressed with the query genes
 * \li retrieve a list of datasets where this coexpression with the query is found to be occurring.
 *
 * \section sec_usage Usage
 * 
 * \subsection ssec_usage_basic Basic Usage
 * 
 * \code
 * SeekServer -t <port> -x <dset_platform_map> -i <gene_map> -d <db_dir> -p <prep_dir> -P <platform_dir>
 * -Q <quant> -n <num_db> -u <sinfo_dir>
 * \endcode
 * 
 * This starts an instance of SeekServer on the indicated port and begins accepting client requests.
 *
 * \subsubsection ssec_cl Client Request Format
 *
 * When a client request comes in, SeekServer looks for the following sequence of 4 strings in the request message:
 *
 * \li \c strSearchDataset. Dataset names, as referred by the \c dset_platform_map, to be used for the search.
 * Delimited by " ".
 *
 * \li \c strQuery. Query gene names, as referred by the \c gene_map, separated by " ".
 *
 * \li \c strOutputDir. Output directory where intermediate results are generated. Must be a directory that the running user of
 * SeekServer has access to. \c /tmp is recommended.
 *
 * \li \c strSearchParameter. A string of the form "1_2_3_4" where each number denotes the following: <br>
 * 1 - the search method, one of \c RBP, \c OrderStatistics, \c EqualWeighting. Recommended \c RBP (also known as the CV weighting). <br>
 * 2 - rbp parameter \a p (a \c float 0.90 - 0.99). Recommended 0.99. <br>
 * 3 - minimum fraction of query required to score each dataset (0 - 1.0). Recommended 0 (no minimum). <br>
 * 4 - distance measure, one of \c Correlation, \c Zscore, \c ZscoreHubbinessCorrected. Recommended \c ZscoreHubbinessCorrected.<br>
 *
 * See Sleipnir::CSeekNetwork for the specification of an incoming string message.
 *
 * Once SeekServer correctly receives the above 4 strings, a search instance using the provided search parameters will
 * be initiated on the server side.
 *
 *
 * \subsubsection ssec_out Outgoing Message Format
 *
 * Each outgoing message is generated upon finishing searching the client's query. In general, if the search is successful,
 * SeekServer will send to the clients these two arrays in sequence:
 * \li a binary \c float array of <b>dataset weights</b>, indicating how datasets are related to the query.
 * \li a binary \c float array of <b>gene scores</b>, indicating how genes are coexpressed with the query. 
 *
 * See Sleipnir::CSeekNetwork for the specification of an outgoing float array.
 *
 *
 * \subsubsection ssec_search Query-independent search setting files and directories
 *
 * These include the following: \c dset_platform_map, \c gene_map, \c db_dir, \c prep_dir, \c platform_dir, \c quant,
 * \c sinfo_dir.
 * For a discussion of these files and directories, please refer to the \ref SeekMiner page in section:
 * Query-independent search setting files and directories.
 *
 *
 * \subsection ssec_usage_detailed Detailed Usage
 * 
 * \include SeekServer/SeekServer.ggo
 * 
 */

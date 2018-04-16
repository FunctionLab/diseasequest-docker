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
 * \page Data2DB Data2DB
 * 
 * Data2DB converts a collection of DAT/DAB files (Sleipnir::CDat) into a simple flatfile database
 * (Sleipnir::CDatabase).  DAT/DAB files organize data so that values for all gene pairs within a single
 * dataset can be accessed efficiently; database files organize data so that values from all datasets for
 * a single gene or gene pair can be accessed efficiently.  This is critical for real-time Bayesian inference (e.g., by \ref BNServer) and for Seek coexpression search
 * (e.g. by \ref SeekMiner, \ref SeekServer).
 * 
 * \section sec_usage Usage
 * 
 * \subsection ssec_usage_basic Basic Usage
 * 
 * \code
 * Data2DB -n <classifier.xdsl> -i <genes.txt> -d <data_dir> -D <database_dir>
 *
 * \endcode
 * Construct a Sleipnir::CDatabase in the directory \c database_dir containing the data from DAT/DAB files
 * in \c data_dir corresponding to nodes in the Bayesian network \c classifier.xdsl and organized using the
 * gene index/name pairs in \c genes.txt (identical in format to \ref Data2Sql).  If many datasets are
 * being processed or the target genome is large, blocking should be used (\c -b and \c -B).
 *
 *
 * \code
 * Data2DB -x <dataset_file_list.txt> -i <gene_map.txt> -D <database_dir> 
 * \endcode
 * Construct a Sleipnir::CDatabase containing the data from DAB files that 
 * are specified in the \c dataset_file_list.txt. The genes are indexed according
 * to \c gene_map.txt. By default, there would be 1000 Sleipnir::CDatabaselet's (DB files)
 * generated, with each containing \a N / 1000 genes. Users can control
 * the number of generated DB files (and indirectly the number of genes contained in each DB)
 * using the \c -f option.
 * 
 * 
 * \subsection ssec_usage_detailed Detailed Usage
 * 
 * \include Data2DB/Data2DB.ggo
 * 
 * <table><tr>
 *	<th>Flag</th>
 *	<th>Default</th>
 *	<th>Type</th>
 *	<th>Description</th>
 * </tr><tr>
 *  <td>-x</td>
 *  <td>None</td>
 *  <td>Dataset file list</td>
 *  <td>A simple one-column listing of path of DAB files. Dataset order in the
 *   CDatabase will correspond to the order in this file. Either this option or the
 *   \c -n option must be specified.</td>
 * </tr><tr>
 *	<td>-n</td>
 *	<td>None</td>
 *	<td>(X)DSL file</td>
 *	<td>Naive Bayesian classifier for which output database will be optimized.  Dataset order in the output
 *		database will correspond to the Bayes net's node order, and the node IDs will be used to load
 *		input DAT/DABs from \c -d. Either this option or the -c -x option must be
 *     specified.</td>
 * </tr><tr>
 *	<td>-i</td>
 *	<td>stdin</td>
 *	<td>Text file</td>
 *	<td>Tab-delimited text file containing two columns, numerical gene IDs (one-based) and unique gene
 *		names (matching those in the input DAT/DAB files).</td>
 * </tr><tr>
 *	<td>-d</td>
 *	<td>.</td>
 *	<td>Directory</td>
 *	<td>Input directory containing DAT/DAB files with names corresponding to the given Bayes net node IDs.</td>
 * </tr><tr>
 *	<td>-D</td>
 *	<td>.</td>
 *	<td>Directory</td>
 *	<td>Output directory in which database files will be stored.</td>
 * </tr><tr>
 *	<td>-f</td>
 *	<td>1000</td>
 *	<td>Integer</td>
 *	<td>Number of separate database files to store in the output directory</td>
 * </tr><tr>
 *	<td>-b</td>
 *	<td>-1</td>
 *	<td>Integer</td>
 *	<td>Number of output files (and hence genes) to process per block.  -1 indicates that all output files
 *		should be created in a single pass.</td>
 * </tr><tr>
 *	<td>-B</td>
 *	<td>-1</td>
 *	<td>Integer</td>
 *	<td>Number of input files (datasets) to process per block.  -1 indicates that all input files should be
 *		read into memory simultaneously.</td>
 * </tr><tr>
 *	<td>-u</td>
 *	<td>off</td>
 *	<td>Flag</td>
 *	<td>If on, buffer each database file in memory during modification and write as a single unit on
 *		completion.  Could in theory speed up database construction on certain disks/filesystems.</td>
 * </tr><tr>
 *	<td>-m</td>
 *	<td>off</td>
 *	<td>Flag</td>
 *	<td>If given, memory map the input files when possible.  DAT and PCL inputs cannot be memmapped.</td>
 * </tr><tr>
 *	<td>-N</td>
 *	<td>off</td>
 *	<td>Flag</td>
 *	<td>If enabled, use Nibble (4 bits) to represent each element rather than the default 8 bits (or a byte).</td>
 * </tr></table>
 */

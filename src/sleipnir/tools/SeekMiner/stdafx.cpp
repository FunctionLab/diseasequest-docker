/*****************************************************************************
* This file is provided under the Creative Commons Attribution 3.0 license.
*
* You are free to share, copy, distribute, transmit, or adapt this work
* PROVIDED THAT you attribute the work to the authors listed below.
* For more information, please see the following web page:
* http://creativecommons.org/licenses/by/3.0/
*
* This file is a part of SEEK (Search-based exploration of expression compendium)
* which is authored and maintained by: Qian Zhu (qzhu@princeton.edu)
*
* If you use this file, please cite the following publication:
* Qian Zhu, Aaron K Wong, Arjun Krishnan, Miriam R Aure, Alicja Tadych, 
* Ran Zhang, David C Corney, Casey S Greene, Lars A Bongo, 
* Vessela N Kristensen, Moses Charikar, Kai Li & Olga G Troyanskaya
* "Targeted exploration and analysis of large cross-platform human 
* transcriptomic compendia" Nat Methods (2015)
*
* This file is a component of the Sleipnir library for functional genomics,
* authored by:
* Curtis Huttenhower (chuttenh@princeton.edu)
* Mark Schroeder
* Maria D. Chikina
* Olga G. Troyanskaya (ogt@princeton.edu, primary contact)
*
* If you use this library for development, or use any other Sleipnir executable
* tools, please also cite the following publication:
* Curtis Huttenhower, Mark Schroeder, Maria D. Chikina, and
* Olga G. Troyanskaya.
* "The Sleipnir library for computational functional genomics"
*****************************************************************************/
#include "stdafx.h"

/*!
 * \page SeekMiner SeekMiner
 *
 * SeekMiner is the main program for integrating coexpressions among thousands of
 * microarray datasets. Users supply the program with a set of genes as input (or query) and the program
 * returns other similar genes with a coexpression to the input genes.
 *
 * The main challenge in performing the user's query is finding the right datasets.
 * As not all microarrays are relevant to exploring the query's coexpression,
 * SeekMiner particularly favors those datasets where the query genes are highly
 * correlated among each other. The intuition is based on the observation that the coregulation
 * between the query genes would suggest that they participate in the same biological process,
 * the biological process involving these genes is highly active. So datasets that pass this criteria would
 * be very informative to the search process.
 *
 * In addition to the default coregulation based weighting, SeekMiner supports other methods
 * of scoring datasets, such as rank-based methods (order statistics) and equal-weighting.
 *
 * Users can easily compare between methods, adjust parameters in the search algorithms, specify
 * the datasets to be integrated, and test a number of different queries with varying length,
 * in order to achieve their desired results.
 *
 * \section sec_usage Usage
 * 
 * \subsection ssec_usage_basic Basic Usage
 * 
 * \code
 * SeekMiner -x <dset_platform_map> -i <gene_map> -q <query> -P <platform_dir> -p <prep_dir> -n <num_db>
 * -d <db_dir> -Q <quant> -o <output_dir> -V <weight_method> -z <distance_measure> -m [-D <search_dset>]
 * \endcode
 * This performs the coexpression search for a list of queries,
 * and outputs the gene-ranking and the dataset weights in the \c output_dir.
 *
 *
 * \subsubsection sec_weight Weighting Datasets
 *
 * SeekMiner supports the following weighting methods (\c -V):
 * \li Query cross-validated weighting (\c CV, default), where we iteratively use a subset of the
 * query to construct a search instance to retrieve the remaining query genes. This is a form of measuring
 * the coregulation of query genes using a cross-validation setup.
 * \li Equal weighting (\c EQUAL), where all datasets are weighted equally.
 * \li Order statistics integration (\c ORDER_STAT), which is outlined in Adler et al (2009).
 * This method computes a P-value statistics by comparing the rank of correlation across datasets to the
 * ranks that would have been generated a null distribution (where correlations are assumed to be
 * randomly scattered and all ranks are equally likely).
 * 
 * The use of \c -V \c CV is highly recommended.
 *
 * \subsubsection sec_distance Distance Measure and Transformations
 *
 * Users can select between Pearson correlations (\c -z \c pearson) or z-scores of Pearson (\c -z \c z_score).
 * Z-scores is the recommended choice because it normalizes the correlation distribution to a standard normal
 * distribution that can be compared across datasets. In addition, SeekMiner provides the following
 * transformations on z-scores to further allow boosting of signals:
 *
 * \li \c --score_cutoff. Cuts off z-scores at a specified value. Z-scores that fall below the cut-off are assigned zero.
 * \li \c --norm_subavg. Subtracts each gene's average z-score. This prevents highly connected genes from being constantly returned with top ranks in the ranking.
 * \li \c --norm_subavg_plat. Normalizes z-score by subtracting the average across the platform and dividing by its standard deviation.
 * This is designed to handle potential platform biases on the z-scores.
 * \li \c --square_z. Squaring the z-score. This is another way to boost the highly correlated gene-pairs.
 * 
 * It is highly recommended to enable \c --norm_subavg.
 *
 * \subsubsection sec_search Search Datasets
 *
 * Users may also define the datasets that they wish to use for integrations in a query-specific way, using \c -D argument.
 * If this argument is absent, all datasets in the compendium will be integrated.
 * If \c -D is used, the search datasets must be selected from the available
 * datasets defined in \c dset_platform_map.
 *
 * \subsubsection sec_output Output
 *
 * The output files are divided according to queries.
 * Starting with the first query (with a file name 0), its final results
 * will consist of three files: \c 0.query, \c 0.dweight, \c 0.gscore.
 * \li The file base name (0) indicates the query index in the list.
 * \li The \c 0.query stores the space-delimited query gene-set in text.
 * \li The \c 0.dweight stores the weightings of datasets as a binary one-dimensional float vector
 * (see \ref SeekEvaluator for displaying a DWEIGHT extension file).
 * \li The \c 0.gscore stores the gene scores as a binary one-dimensional float vector
 * (see \ref SeekEvaluator for displaying a GSCORE extension file).
 *
 * \subsubsection sec_files Query-independent search setting files and directories
 *
 * \c -x \c dset_platform_map
 *
 * Tab-delimited text file containing two columns, the dataset name,
 * and the corresponding platform name. Below is a few sample lines:
 * \code
 * GSE15913.GPL570.pcl  GPL570
 * GSE16122.GPL2005.pcl GPL2005
 * GSE16797.GPL570.pcl  GPL570
 * GSE16836.GPL570.pcl  GPL570
 * GSE17351.GPL570.pcl  GPL570
 * GSE17537.GPL570.pcl  GPL570
 * \endcode
 * Note that although the dataset name looks like a file name, it does not
 * need to be a valid file name, as long as it properly and uniquely describes
 * the dataset. Here, the dataset is uniquely identified by a GSE ID and a GPL ID
 * combination. In addition, the ordering of the datasets in this file must match
 * the order of the datasets in the CDatabaselet (ie DB files).
 *
 * \c -i \c gene_map
 *
 * Tab-delimited gene-map file. Maps the genes to an ID between 0 to N where N is
 * the genome size. Example:
 * \code
 * 1    1
 * 2    10
 * 3    100
 * 4    1000
 * 5    10000
 * 6    100008589
 * 7    100009676
 * 8    10001
 * 9    10002
 * 10   10003
 * 11   100033413
 * 12   100033414
 * \endcode
 * The ordering of the genes in this file must match the order of genes
 * in the CDatabaselets (DB files).
 *
 * \c -q \c query
 *
 * The file can contain multiple queries that are listed one query per line.
 * The genes in each query are separated by spaces. Example:
 * \code
 * 10003 10002 10001
 * 634 6265
 * \endcode
 * The names of the genes must be selected from the genes in the \c gene_map.
 * The maximum length of the query depends on the amount of available memory in the system.
 * It is recommended to keep each query less than 100 genes.
 *
 * \c -D \c search_dset
 *
 * This file defines the list of datasets to be used for the query coexpression search.
 * The file is defined in a query specific way.
 * An example is provided below:
 * \code
 * GSE15913.GPL570.pcl GSE16122.GPL2005.pcl GSE16836.GPL570.pcl ...
 * GSE14933.GPL570.pcl GSE15162.GPL2005.pcl GSE15566.GPL570.pcl ...
 * ...
 * \endcode
 * where each line, corresponding to a query, is a space-separated dataset list for the query.
 * The dataset names must be selected from the file \c dset_platform_map.
 *
 * \c -P \c platform_dir
 *
 * Directory that contains the following 3 files:
 * \li \c all_platforms.gplatavg. the platform average z-scores
 * \li \c all_platforms.gplatstdev. the platform z-score standard deviation
 * \li \c all_platforms.gplatorder. the order of platforms
 *
 * These binary files are generated by \ref SeekPrep. The specification of this directory is
 * necessary for \c --norm_subavg_plat.
 *
 * \c -p \c prep_dir
 *
 * Directory that contains the gene presence files and the gene average files:
 * \li Gene presence (GPRES files): indicates the presence/absence of genes in a dataset
 * \li Gene average (GAVG files): indicates the average z-score of each gene in a dataset
 *
 * There should be one pair of these files for <b>every</b> dataset that is specified
 * in \c dset_platform_map. Generated by \ref SeekPrep.
 *
 * \c -d \c db_dir
 *
 * Directory that contains the CDatabase (all of the DB files).
 *
 * \c -Q \c quant
 *
 * The \c quant file specifies how the z-scores are binned. This is necessary for properly reading
 * the z-scores, because the z-scores are stored as binned values on disk. This quant file is used
 * to convert them back to z-scores when they are read from disk.
 * Currently, the maximum number of bins supported is 255.
 * A snapshot of the \c quant file is below:
 * \code
 * -5.00 -4.96 -4.92 -4.88 -4.84 -4.80 -4.76 -4.72 -4.68 -4.64 -4.60 -4.56 -4.52 ...
 * \endcode
 * The bin boundaries are separated by spaces.
 *
 * \c -o \c output_dir
 *
 * Directory that will contain the search results.
 * 
 * \c -u \c sinfo_dir
 *
 * Directory that contains the SINFO files, which list a dataset's average z-score between all pairs of genes
 * and the standard deviation. If this directory is provided, there should be one SINFO file for <b>
 * every</b> dataset in \c dset_platform_map. Generated by \ref SeekPrep.
 *
 * 
 * \subsection ssec_usage_detailed Detailed Usage
 * 
 * \include SeekMiner/SeekMiner.ggo
 * 
 */

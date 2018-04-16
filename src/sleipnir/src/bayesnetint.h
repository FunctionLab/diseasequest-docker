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
#ifndef BAYESNETINT_H
#define BAYESNETINT_H

#include "fullmatrix.h"

namespace Sleipnir {

class CDat;
class CPCL;
class CPCLPair;
class IDataset;

/*!
 * \brief
 * Encapsulates a Bayesian network with arbitrary structure and node types.
 * 
 * IBayesNet provides an interface for Bayesian graphical models.  These can have an arbitrary graph
 * structure, and implementations of the interface can provide arbitrary node types.  Inference and
 * parameter learning functions are exposed that operate directly on Sleipnir datatypes such as IDataset,
 * CDat, and CPCL.  The only strict requirement of nodes is that they provide unique string labels and
 * can expose their parameters by way of a CDataMatrix, although the semantics of those parameters are
 * not constrained.
 */
class IBayesNet {
public:
	/*!
	 * \brief
	 * Load a Bayes net from a file.
	 * 
	 * \param szFile
	 * Path to file.
	 * 
	 * \returns
	 * True if Bayes net was loaded succesfully.
	 * 
	 * \remarks
	 * Specific behavior is implementation specific; it is assumed that the network will be completely
	 * reinitialized from the given file, although it may be left in an inconsistent state if the return
	 * value is false.
	 */
	virtual bool Open( const char* szFile ) = 0;
	/*!
	 * \brief
	 * Save a Bayes net to a file.
	 * 
	 * \param szFile
	 * Path to file.
	 * 
	 * \returns
	 * True if Bayes net was saved succesfully.
	 * 
	 * \remarks
	 * Specific behavior is implementation specific; the Bayes net will not be modified, but the contents
	 * of the output file may be inconsistent if the return value is false.
	 */
	virtual bool Save( const char* szFile ) const = 0;
	/*!
	 * \brief
	 * Learn conditional probabilities from data using Expectation Maximization, naive Bayesian learning, or
	 * Extended Logistic Regression.
	 * 
	 * \param pDataset
	 * Dataset to be used for learning.
	 * 
	 * \param iIterations
	 * Maximum number of iterations for EM or ELR.
	 * 
	 * \param fZero
	 * If true, assume all missing values are zero (i.e. the first bin).
	 * 
	 * \param fELR
	 * If true, use ELR to learn network parameters.
	 * 
	 * \returns
	 * True if parameters were learned successfully.
	 * 
	 * Using the given IDataset, learn parameters for the underlying Bayes network.  If requested, learning
	 * is performed discriminatively using Extended Logistic Regression due to Greiner, Zhou, et al.
	 * Otherwise, maximum likelihood estimates are used for naive structures, and Expectation Maximization
	 * is used for other network structures.
	 * 
	 * \remarks
	 * The order of datasets in the given IDataset must correspond to the order of nodes within the Bayes
	 * network, and the first dataset (index 0) is assumed to be a gold standard.  Only data for which
	 * IDataset::IsExample is true will be used, which usually means that the first dataset and at least
	 * one other dataset must have a value.
	 */
	virtual bool Learn( const IDataset* pDataset, size_t iIterations, bool fZero = false,
		bool fELR = false ) = 0;
	/*!
	 * \brief
	 * Perform Bayesian inference to obtain probabilities for each element of a dataset.
	 * 
	 * \param pDataset
	 * Dataset to be used as input for inference.
	 * 
	 * \param vecvecdResults
	 * Vector of output probabilities; each element of the outer vector represents the result for one
	 * gene pair, and each element of the inner vectors represents the probability for one possible value
	 * from the output node (i.e. the answer).
	 * 
	 * \param fZero
	 * If true, assume all missing values are zero (i.e. the first bin).
	 * 
	 * \returns
	 * True if evaluation was successful.
	 * 
	 * The inverse of the corresponding IBayesNet::Learn method; given an IDataset, ignore the first
	 * (gold standard) dataset and infer the corresponding output probabilities for each other gene pair
	 * for which data is available.  For each gene pair within the IDataset for which IDataset::IsExample is
	 * true, vecvecdResults will contain one vector.  This vector will contain inferred probabilities for
	 * each possible value of the output node, generally the probability of functional unrelatedness (i.e.
	 * one minus the probability of functional relationship).
	 * 
	 * \remarks
	 * The order of datasets in the given IDataset must correspond to the order of nodes within the Bayes
	 * network, and the first dataset (index 0) is assumed to be a gold standard (and is thus ignored).  Only
	 * data for which IDataset::IsExample is true will be used, which usually means that at least one other
	 * dataset must have a value.  If the output node can take N values, each output vector will contain only
	 * the first N-1 probabilities, since the Nth can be calculated to sum to one.
	 */
	virtual bool Evaluate( const IDataset* pDataset, std::vector<std::vector<float> >& vecvecdResults,
		bool fZero = false ) const = 0;
	/*!
	 * \brief
	 * Perform Bayesian inference to obtain probabilities for each element of a dataset.
	 * 
	 * \param pDataset
	 * Dataset to be used as input for inference.
	 * 
	 * \param DatResults
	 * Description of parameter DatResults.
	 * 
	 * \param fZero
	 * If true, assume all missing values are zero (i.e. the first bin).
	 * 
	 * \returns
	 * True if evaluation was successful.
	 * 
	 * The inverse of the corresponding IBayesNet::Learn method; given an IDataset, ignore the first
	 * (gold standard) dataset and infer the corresponding output probability for each other gene pair
	 * for which data is available.  For each gene pair within the IDataset for which IDataset::IsExample is
	 * true, the probability of functional relationship (i.e. the largest possible value of the output node)
	 * will be placed in the given CDat.
	 * 
	 * \remarks
	 * The order of datasets in the given IDataset must correspond to the order of nodes within the Bayes
	 * network, and the first dataset (index 0) is assumed to be a gold standard (and is thus ignored).  Only
	 * data for which IDataset::IsExample is true will be used, which usually means that at least one other
	 * dataset must have a value.
	 */
	virtual bool Evaluate( const IDataset* pDataset, CDat& DatResults, bool fZero = false ) const = 0;
	/*!
	 * \brief
	 * Perform Bayesian inference to obtain probabilities given values for each other Bayes net node.
	 * 
	 * \param vecbDatum
	 * One-indexed values for each node in the Bayes net (zero indicates missing data).
	 * 
	 * \param vecdResults
	 * Inferred probabilities for each possible value of the requested node.
	 * 
	 * \param fZero
	 * If true, assume all missing values are zero (i.e. the first bin).
	 * 
	 * \param iNode
	 * The node for which output probabilities are inferred.
	 * 
	 * \param fIgnoreMissing
	 * If true, do not default missing values to zero or any other value.
	 * 
	 * \returns
	 * True if evaluation was successful.
	 * 
	 * This Evaluate assumes a discrete Bayes net and, given a vector of evidence values for each node,
	 * infers the probability distribution over possible values of the requested node.  Note that vecbDatum
	 * contains <b>one plus</b> the discrete bin value of each node, and a value of zero indicates missing
	 * data for the corresponding node.
	 * 
	 * \remarks
	 * vecbDatum should contain one plus the discrete bin value of each node, and a value of zero indicates
	 * missing data for the corresponding node.  If the requested output node can take N values, the output
	 * vector will contain only the first N-1 probabilities, since the Nth can be calculated to sum to one.
	 */
	virtual bool Evaluate( const std::vector<unsigned char>& vecbDatum, std::vector<float>& vecdResults,
		bool fZero = false, size_t iNode = 0, bool fIgnoreMissing = false ) const = 0;
	/*!
	 * \brief
	 * Perform Bayesian inference to obtain probabilities over all nodes in the network given some amount of data.
	 * 
	 * \param PCLData
	 * Input data; each column (experiment) is mapped by label to a node in the Bayes net, and PCL entries
	 * correspond to observed (or missing) data values.
	 * 
	 * \param PCLResults
	 * Output probabilities; each column (experiment) is mapped to a node:value pair from the Bayes net, and
	 * PCL entries correspond to the probability of that value in that node.
	 * 
	 * \param fZero
	 * If true, assume all missing values are zero (i.e. the first bin).
	 * 
	 * \param iAlgorithm
	 * Implementation-specific ID of the Bayesian inference algorithm to use.
	 * 
	 * \returns
	 * True if evaluation was successful.
	 * 
	 * This version of Evaluate will perform one Bayesian inference for each row (gene) of the given PCLData.
	 * Here, each PCL "experiment" column corresponds to a node in the Bayes net as identified by the
	 * experiment labels in the PCL and the IDs of the Bayes net nodes.  Values are read from the given PCL
	 * and (if present; missing values are allowed) discretized into Bayes net value bins using the
	 * accompanying quantization information.  For each input row, all given non-missing values are observed
	 * for the appropriate Bayes net nodes, and Bayesian inference is used to provide probabilities for
	 * each remaining, unobserved node value.
	 * 
	 * \remarks
	 * PCLResults must be initialized with the correct number of experimental columns before calling
	 * Evaluate; that is, the total number of node values in the Bayes net.  For example, if the Bayes net
	 * has three nodes A, B, and C, node A can take two values 0 and 1, and nodes B and C can take values
	 * 0, 1, and 2, then PCLResults must have 8 experimental columns corresponding to A:0, A:1, B:0, B:1, B:2,
	 * C:0, C:1, and C:2.  Columns of PCLData are mapped to Bayes net nodes by experiment and node labels;
	 * experiment labels not corresponding to any Bayes net node ID are ignored, and Bayes net nodes with no
	 * corresponding experiment are assumed to be unobserved (hidden).  Only the genes in PCLResults are used,
	 * and they need not be in the same order as in PCLData.
	 */
	virtual bool Evaluate( const CPCLPair& PCLData, CPCL& PCLResults, bool fZero = false,
		int iAlgorithm = -1 ) const = 0;
	/*!
	 * \brief
	 * Retrieve the string IDs of all nodes in the Bayes net.
	 * 
	 * \param vecstrNodes
	 * Output containing the IDs of all nodes in the Bayes net.
	 */
	virtual void GetNodes( std::vector<std::string>& vecstrNodes ) const = 0;
	/*!
	 * \brief
	 * Returns the number of different values taken by the requested node.
	 * 
	 * \param iNode
	 * Bayes net node for which values should be returned.
	 * 
	 * \returns
	 * Number of different values taken by the requested node.
	 * 
	 * \remarks
	 * Not applicable for continuous nodes.
	 */
	virtual unsigned char GetValues( size_t iNode ) const = 0;
	/*!
	 * \brief
	 * Returns true if any node in the Bayes net is non-discrete (e.g. Gaussian, etc.)
	 * 
	 * \returns
	 * True if any node in the Bayes net is continuous.
	 */
	virtual bool IsContinuous( ) const = 0;
	/*!
	 * \brief
	 * Returns true if the requested node is non-discrete (e.g. Gaussian, etc.)
	 * 
	 * \param iNode
	 * Node to be inspected.
	 * 
	 * \returns
	 * True if the requested node is continuous.
	 */
	virtual bool IsContinuous( size_t iNode ) const = 0;
	/*!
	 * \brief
	 * Randomizes every parameter in the Bayes net.
	 * 
	 * \remarks
	 * Parameter values are generated uniformly at random and normalized to represent a valid probability
	 * distribution.
	 */
	virtual void Randomize( ) = 0;
	/*!
	 * \brief
	 * Randomizes every parameter the requested node.
	 * 
	 * \param iNode
	 * Index of node to be randomized.
	 * 
	 * \remarks
	 * Parameter values are generated uniformly at random and normalized to represent a valid probability
	 * distribution.
	 */
	virtual void Randomize( size_t iNode ) = 0;
	/*!
	 * \brief
	 * Reverses the parameters of the requested node over its possible values.
	 * 
	 * \param iNode
	 * Index of node to be reversed.
	 * 
	 * "Vertically" reverses the parameters of the requested node.  That is, if the requested node can take
	 * values 0 through 3, then for each setting of the parents' values, Pnew(0|parents) = Pold(3|parents),
	 * Pnew(1|parents) = Pold(2|parents), Pnew(2|parents) = Pold(1|parents), and Pnew(3|parents) =
	 * Pold(0|parents).
	 * 
	 * \remarks
	 * May be ignored by some implementations, particularly continuously valued nodes.
	 */
	virtual void Reverse( size_t iNode ) = 0;
	/*!
	 * \brief
	 * Retrieves the parameters of the requested Bayes net node.
	 * 
	 * \param iNode
	 * Index of node for which parameters should be retrieved.
	 * 
	 * \param MatCPT
	 * Parameters of the requested node in tabular form; the columns of the matrix represent parental
	 * values, the rows node values.
	 * 
	 * \returns
	 * True if parameter retrieval succeeded, false if it failed or the requested node has more than one
	 * parent.
	 * 
	 * Retrieves node parameters in an implementation-specific manner, often only allowing nodes with
	 * at most one parent.  For discrete nodes, matrix entries are generally conditional probabilities.
	 * For continuous nodes, matrix entries may represent distribution parameters such as Gaussian mean
	 * and standard deviation.
	 * 
	 * \remarks
	 * Only allowed for nodes with at most one parent; nodes with more parents are supported by some
	 * implementations, but their parameters can't be retrieved by this function.
	 */
	virtual bool GetCPT( size_t iNode, CDataMatrix& MatCPT ) const = 0;
};

}

#endif // BAYESNETINT_H

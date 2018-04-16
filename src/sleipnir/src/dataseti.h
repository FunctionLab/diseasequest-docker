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
#ifndef DATASETI_H
#define DATASETI_H

#include <set>

#include "examplei.h"
#include "fullmatrix.h"
#include "halfmatrix.h"

namespace Sleipnir {

class CCompactMatrix;
class CDataPair;
class IBayesNet;
class IDataset;

class CDataImpl {
	friend class CDataFilter;
	friend class CDataMask;
protected:
	static const char	c_cSeparator	= '/';
	static const char	c_szDat[];
	static const char	c_szDab[];

	static void FilterGenes( IDataset*, const CGenes&, CDat::EFilter );

	size_t OpenMax( const char*, const std::vector<std::string>&, bool,
		std::vector<std::string>&, std::set<std::string>* = NULL );
	bool OpenGenes( std::istream&, bool, bool, std::set<std::string>& ) const;
	bool OpenGenes( const std::vector<std::string>& );

	size_t GetGene( const std::string& ) const;
	bool OpenBinary( std::istream& );
	const unsigned char* OpenBinary( const unsigned char* );
	void SaveBinary( std::ostream& ) const;

	bool IsHidden( size_t iNode ) const {

		return ( m_veciMapping[ iNode ] == -1 ); }

	const std::string& GetGene( size_t iGene ) const {

		return m_vecstrGenes[ iGene ]; }

	size_t GetGenes( ) const {

		return m_vecstrGenes.size( ); }

	const std::vector<std::string>& GetGeneNames( ) const {

		return m_vecstrGenes; }

	size_t GetExperiments( ) const {

		return m_veciMapping.size( ); }

	unsigned char GetBins( size_t iExp ) const {

		return m_veccQuants[ iExp ]; }

	bool						m_fContinuous;
	std::vector<size_t>			m_veciMapping;
	std::vector<std::string>	m_vecstrGenes;
	std::vector<unsigned char>	m_veccQuants;
};

class CDatasetImpl : protected CDataImpl {
protected:
	CDatasetImpl( );
	virtual ~CDatasetImpl( );

	void Reset( );
	bool Open( const CDataPair&, size_t );
	bool Open( const CDataPair*, const char*, const IBayesNet* );
	void SaveText( std::ostream& ) const;
	void SaveBinary( std::ostream& ) const;
	float GetContinuous( size_t, size_t, size_t ) const;

	void**	m_apData;
};

class CDataOverlayImpl {
protected:
	CDataOverlayImpl( ) : m_pDataset(NULL) { }

	const std::vector<std::string>& GetGeneNames( ) const;
	size_t GetExperiments( ) const;
	size_t GetGene( const std::string& ) const;
	size_t GetBins( size_t ) const;
	size_t GetGenes( ) const;
	bool IsHidden( size_t ) const;
	size_t GetDiscrete( size_t, size_t, size_t ) const;
	float GetContinuous( size_t, size_t, size_t ) const;
	const std::string& GetGene( size_t ) const;
	void Save( std::ostream&, bool ) const;

	const IDataset*	m_pDataset;
};

class CDataMaskImpl : protected CDataOverlayImpl {
protected:
	CBinaryMatrix	m_Mask;
};

class CDataFilterImpl : protected CDataOverlayImpl {
protected:
	CDataFilterImpl( ) : m_pGenes(NULL), m_pAnswers(NULL) { }

	const CGenes*		m_pGenes;
	CDat::EFilter		m_eFilter;
	std::vector<bool>	m_vecfGenes;
	const CDat*			m_pAnswers;
	std::vector<size_t>	m_veciAnswers;
};

class CDataSubsetImpl : protected CDataImpl {
protected:
	bool Open( const CDataPair&, size_t );

	size_t						m_iSize;
	size_t						m_iOffset;
	std::vector<std::string>	m_vecstrData;
	CFullMatrix<CExampleImpl>	m_Examples;
};

class CDatasetCompactImpl : protected CDataImpl {
protected:
	CDatasetCompactImpl( );
	virtual ~CDatasetCompactImpl( );

	bool Open( const CDataPair&, size_t );
	bool Open( const char*, const IBayesNet*, const CGenes* = NULL, const CGenes* = NULL );
	bool Open( const unsigned char* );
	virtual void Remove( size_t, size_t );
	size_t GetDiscrete( size_t, size_t, size_t ) const;
	void SaveText( std::ostream& ) const;
	void SaveBinary( std::ostream& ) const;
	virtual bool IsExample( size_t, size_t ) const;

	uint32_t		m_iData;
	CCompactMatrix*	m_aData;
};

}

#endif // DATASETI_H

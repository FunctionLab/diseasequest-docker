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
#ifndef DATI_H
#define DATI_H

#include <map>

#include "halfmatrix.h"
#include "file.h"
#include "measure.h"
#include "meta.h"
#include "pcl.h"

namespace Sleipnir {

class CColor;
class CGenes;
class CSlim;

class CDatImpl : protected CFile {
protected:
	typedef std::map<std::string,size_t>	TMapStrI;
	typedef std::vector<std::string>		TVecStr;
	typedef std::vector<float>				TAF;
	typedef std::vector<TAF>				TAAF;

	static const size_t		c_iGeneLimit	= 100000;
	static const size_t		c_iNeighborhood	= 40;
	static const size_t		c_iDegree		= 1;
	static const char		c_acComment[];
	static const CColor&	c_ColorMid;
	static const CColor&	c_ColorMin;
	static const CColor&	c_ColorMax;

	static size_t MapGene( TMapStrI&, TVecStr&, const std::string& );
	static void ResizeNaN( TAF&, size_t );
	static void DabGene( std::istream&, char* );

	CDatImpl( ) : m_pMeasure(NULL), m_pPCL(NULL), m_abData(NULL), m_iData(0), m_aadData(NULL), m_hndlData(0) { }
	virtual ~CDatImpl( );

	void Reset( );
	bool OpenPCL( std::istream&, size_t, bool );
	bool OpenText( std::istream&, float, bool );
	bool OpenBinary( std::istream&, bool = false );
	bool OpenSparse( std::istream& );
	bool OpenQdab( std::istream& );
	bool OpenGenes( std::istream&, bool, bool );
	void SaveText( std::ostream& ) const;
	void SaveBinary( std::ostream& ) const;
	void SaveSparse( std::ostream& ) const;
	void SaveGenes( std::ostream& ) const;
	size_t GetGene( const std::string& ) const;
	void SlimCache( const CSlim&, std::vector<std::vector<size_t> >& ) const;
	void AveStd( double&, double&, size_t&, size_t = -1 ) const;
	void NormalizeMinmax( );
	void NormalizeMinmaxNPone( );
	void NormalizePCC( );
	void NormalizeStdev( );
	void NormalizeSigmoid( );
	void NormalizeNormCDF( );
	void OpenHelper( const CGenes*, float );
	void OpenHelper( const CGenes*, const CGenes*, float );
	bool OpenHelper( );
	bool OpenMemmap( const unsigned char* );
	void FilterGenesGraph( const CGenes&, std::vector<bool>&, size_t, float, bool, bool, const std::vector<float>* );

	struct size_t_comp {
	    bool operator ()(size_t const& a, size_t const& b) const {
	    	return (a<b);
	    }
	} size_t_comp;


	float* GetFullRow(const size_t &iY){
		float *d_array = m_Data.GetFullRow(iY);
		d_array[iY] = CMeta::GetNaN();
		return d_array;
	}

	float& Get( size_t iX, size_t iY ) const {
		static float	s_dRet;

		return ( m_pMeasure ? ( s_dRet = (float)m_pMeasure->Measure( m_pPCL->Get( iX ), m_pPCL->GetExperiments( ),
			m_pPCL->Get( iY ), m_pPCL->GetExperiments( ) ) ) :
			( ( iX == iY ) ? ( s_dRet = CMeta::GetNaN( ) ) : m_Data.Get( iX, iY ) ) ); }

	bool Set( size_t iX, size_t iY, float dValue ) {

		if( iX == iY )
			return false;

		m_Data.Set( iX, iY, dValue );
		return true; }

	bool Set( size_t iX, const float* adValues ) {

		m_Data.Set( iX, adValues );
		return true; }

	size_t GetGenes( ) const {

		return ( m_pPCL ? m_pPCL->GetGenes( ) : m_vecstrGenes.size( ) ); }

	size_t GetGeneIndex(const std::string &strGene) const {
		std::map<std::string, size_t>::const_iterator	iterGene;
		return ( ( ( iterGene = m_mapstrGenes.find( strGene ) ) == m_mapstrGenes.end( ) ) ? -1 :
			iterGene->second );
	}

	std::string GetGene( size_t iGene ) const {

		return ( m_pPCL ? m_pPCL->GetGene( iGene ) : m_vecstrGenes[ iGene ] ); }

	const std::vector<std::string>& GetGeneNames( ) const {

		return ( m_pMeasure ? m_pPCL->GetGeneNames( ) : m_vecstrGenes ); }

	void EstimateSeekPositions(istream &istm){
		m_iHeader = istm.tellg();
		size_t i;
		m_veciSeekPos.resize(m_vecstrGenes.size());
		m_veciSeekPos[0] = 0;
		for(i=1; i<m_vecstrGenes.size()-1; i++){
			m_veciSeekPos[i] = m_veciSeekPos[i-1] +
				(sizeof(float)*(m_vecstrGenes.size()-1 - i));
		}
	}

	float* GetRowSeek(std::istream& istm, const std::string &strGene) const;
	float* GetRowSeek(std::istream& istm, const size_t &ind) const;
	bool OpenHeader(std::istream& istm);


	CDistanceMatrix	m_Data;
	TVecStr			m_vecstrGenes;
	std::map<std::string, size_t> m_mapstrGenes;
// PCL back end
	CPCL*			m_pPCL;
	bool			m_fPCLMemory;
	const IMeasure*	m_pMeasure;
	bool			m_fMeasureMemory;
// Memory mapped back end
	unsigned char*	m_abData;
	size_t			m_iData;
	HANDLE			m_hndlData;
	float**			m_aadData;
// Seek positions
	std::vector<size_t>	m_veciSeekPos;
	size_t			m_iHeader;
	bool			m_fSeek;
	/* handle used to open this file
	 * used for reading sparse number of values
	 * without reading the entire file
	 */
	ifstream	m_ifsm;

};

}

#endif // DATI_H

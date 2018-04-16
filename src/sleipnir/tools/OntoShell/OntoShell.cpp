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
#include "cmdline.h"
#include "parserconsole.h"

static const char	c_szCat[]		= "cat";
static const char	c_szFind[]		= "find";
static const char	c_szPrompt[]	= "> ";

static const CParserConsole*	g_pParser;
static bool						g_fGenes;
static bool						g_fZeroes;

bool ProcessLine( const char* );
char** CompletionAll( const char*, int, int );
char* CompletionCommands( const char*, int );
char* CompletionMembers( const char*, int );
size_t CompletionGetParents( const CParser::SLocation&, vector<string>& );
size_t CompletionGetKids( const CParser::SLocation&, vector<string>& );
size_t CompletionGetGenes( const CParser::SLocation&, vector<string>& );

int main( int iArgs, char** aszArgs ) {
	COntologyKEGG			KEGG;
	COntologyOBO				GOBP, GOMF, GOCC;
	COntologyMIPS			MIPS;
	COntologyMIPSPhenotypes	MIPSPhen;
	CGenome					Genome;
	ifstream				ifsmOnto, ifsmGenes;
	gengetopt_args_info		sArgs;
	char*					szLine;
	int						iRet;
	const IOntology*		apOntologies[]
		= { &KEGG, &GOBP, &GOMF, &GOCC, &MIPS, &MIPSPhen, NULL };
	CParserConsole			Parser( apOntologies, Genome );

	g_pParser = &Parser;
	iRet = cmdline_parser2( iArgs, aszArgs, &sArgs, 0, 1, 0 );
	if( sArgs.config_arg )
		iRet = cmdline_parser_configfile( sArgs.config_arg, &sArgs, 0, 0, 1 ) && iRet;
	if( iRet ) {
		cmdline_parser_print_help( );
		return iRet; }
	g_fZeroes = !!sArgs.zeroes_flag;
	CMeta Meta( sArgs.verbosity_arg );

	if( sArgs.features_arg ) {
		ifsmGenes.open( sArgs.features_arg );
		if( !Genome.Open( ifsmGenes ) ) {
			cerr << "Could not open: " << sArgs.features_arg << endl;
			return 1; }
		ifsmGenes.close( ); }

	if( sArgs.kegg_arg ) {
		ifsmOnto.open( sArgs.kegg_arg );
		if( !KEGG.Open( ifsmOnto, Genome, sArgs.kegg_org_arg, !!sArgs.altids_flag ) ) {
			cerr << "Could not open: " << sArgs.kegg_arg << endl;
			return 1; }
		ifsmOnto.close( ); }

	if( sArgs.go_onto_arg ) {
		ifsmOnto.clear( );
		ifsmOnto.open( sArgs.go_onto_arg );
		if( sArgs.go_anno_arg ) {
			ifsmGenes.clear( );
			ifsmGenes.open( sArgs.go_anno_arg ); }
		if( !COntologyOBO::Open( ifsmOnto, ifsmGenes, Genome, GOBP, GOMF, GOCC, !!sArgs.dbids_flag,
			!!sArgs.altids_flag ) ) {
			cerr << "Could not open: " << sArgs.go_onto_arg << ", " << sArgs.go_anno_arg << endl;
			return 1; }
		ifsmOnto.close( );
		if( sArgs.go_anno_arg )
			ifsmGenes.close( ); }

	if( sArgs.mips_onto_arg ) {
		ifsmOnto.clear( );
		ifsmOnto.open( sArgs.mips_onto_arg );
		if( sArgs.mips_anno_arg ) {
			ifsmGenes.clear( );
			ifsmGenes.open( sArgs.mips_anno_arg ); }
		if( !MIPS.Open( ifsmOnto, ifsmGenes, Genome ) ) {
			cerr << "Could not open: " << sArgs.mips_onto_arg << ", " << sArgs.mips_anno_arg <<
				endl;
			return 1; }
		ifsmOnto.close( );
		if( sArgs.mips_anno_arg )
			ifsmGenes.close( ); }

	if( sArgs.mipsp_onto_arg ) {
		ifsmOnto.clear( );
		ifsmOnto.open( sArgs.mipsp_onto_arg );
		if( sArgs.mipsp_anno_arg ) {
			ifsmGenes.clear( );
			ifsmGenes.open( sArgs.mipsp_anno_arg ); }
		if( !MIPSPhen.Open( ifsmOnto, ifsmGenes, Genome ) ) {
			cerr << "Could not open: " << sArgs.mipsp_onto_arg << ", " <<
				sArgs.mipsp_anno_arg << endl;
			return 1; }
		ifsmOnto.close( );
		if( sArgs.mipsp_anno_arg )
			ifsmGenes.close( ); }

	if( sArgs.exec_arg )
		Parser.ProcessLine( sArgs.exec_arg );
	else {
		rl_attempted_completion_function = CompletionAll;
		do {
			if( !( szLine = readline( ( Parser.GetLocation( ).ToString( false ) +
				"> " ).c_str( ) ) ) )
				break;
			add_history( szLine );
			Parser.ProcessLine( szLine );
			free( szLine ); }
		while( true ); }

	return 0; }

char** CompletionAll( const char* szText, int iStart, int iEnd ) {

	if( !iStart )
		return rl_completion_matches( szText, CompletionCommands );
	if( !strncmp( rl_line_buffer, c_szFind, strlen( c_szFind ) ) )
		return NULL; // rl_completion_matches( szText, rl_filename_completion_function );

	g_fGenes = !strncmp( rl_line_buffer, c_szCat, strlen( c_szCat ) );
	return rl_completion_matches( szText, CompletionMembers ); }

char* CompletionCommands( const char* szText, int iState ) {
	static size_t	iTry, iLen;
	const char*	szCur;

	if( !iState ) {
		rl_completion_append_character = ' ';
		iTry = 0;
		iLen = strlen( szText ); }

	while( true ) {
		if( !( szCur = CParser::GetCommand( iTry++ ) ) )
			break;
		if( !strncmp( szCur, szText, iLen ) )
			return _strdup( szCur ); }

	return NULL; }

char* CompletionMembers( const char* szText, int iState ) {
	static size_t				iTry, iLen, iLinks;
	static CParser::SLocation	sLoc;
	static const char*			szReal;
	static vector<string>		vecstrTries;
	static string				strBase;

	if( !iState ) {
		iTry = 0;
		vecstrTries.clear( );
		sLoc = g_pParser->GetLocation( szText, false );
		iLinks = CompletionGetParents( sLoc, vecstrTries );
		iLinks += CompletionGetKids( sLoc, vecstrTries );
		if( g_fGenes )
			CompletionGetGenes( sLoc, vecstrTries );
		strBase = szText;
		if( szReal = strrchr( szText, CParser::c_cSep ) )
			szReal++;
		else
			szReal = szText;
		strBase.resize( szReal - szText );
		iLen = strlen( szReal ); }
	if( !sLoc.IsValid( ) )
		return NULL;

	while( iTry < vecstrTries.size( ) )
		if( !strncmp( szReal, vecstrTries[ iTry++ ].c_str( ), iLen ) ) {
			rl_completion_append_character = ( iTry <= iLinks ) ? '/' : ' ';
			return _strdup( ( strBase + vecstrTries[ iTry - 1 ] ).c_str( ) ); }

	return NULL; }

size_t CompletionGetParents( const CParser::SLocation& sLoc,
	vector<string>& vecstrParents ) {
	const IOntology*	pOnto;
	size_t				i, iSize;

	if( !( pOnto = sLoc.m_pOnto ) || ( sLoc.m_iNode == -1 ) )
		return 0;

	iSize = vecstrParents.size( );
	for( i = 0; i < pOnto->GetParents( sLoc.m_iNode ); ++i )
		vecstrParents.push_back( pOnto->GetID( pOnto->GetParent( sLoc.m_iNode, i ) ) );

	return ( vecstrParents.size( ) - iSize ); }

size_t CompletionGetKids( const CParser::SLocation& sLoc, vector<string>& vecstrKids ) {
	size_t				i, iSize, iChild;
	const IOntology*	pOnto;

	iSize = vecstrKids.size( );
	if( sLoc.m_pOnto && ( sLoc.m_iNode == -1 ) ) {
		for( i = 0; i < sLoc.m_pOnto->GetNodes( ); ++i )
			if( !sLoc.m_pOnto->GetParents( i ) && sLoc.m_pOnto->GetGenes( i, true ) )
				vecstrKids.push_back( sLoc.m_pOnto->GetID( i ) ); }
	else if( sLoc.m_pOnto ) {
		for( i = 0; i < sLoc.m_pOnto->GetChildren( sLoc.m_iNode ); ++i )
			if( sLoc.m_pOnto->GetGenes( iChild = sLoc.m_pOnto->GetChild( sLoc.m_iNode,
				i ), true ) )
				vecstrKids.push_back( sLoc.m_pOnto->GetID( iChild ) ); }
	else
		for( i = 0; i < g_pParser->GetOntologies( ); ++i )
			if( g_pParser->GetGenome( ).CountGenes( pOnto = g_pParser->GetOntology( i ) ) )
				vecstrKids.push_back( pOnto->GetID( ) );

	return ( vecstrKids.size( ) - iSize ); }

size_t CompletionGetGenes( const CParser::SLocation& sLoc, vector<string>& vecstrGenes ) {
	size_t	i, iSize;

	if( !sLoc.m_pOnto || ( sLoc.m_iNode == -1 ) )
		return 0;

	iSize = vecstrGenes.size( );
	for( i = 0; i < sLoc.m_pOnto->GetGenes( sLoc.m_iNode ); ++i )
		vecstrGenes.push_back( sLoc.m_pOnto->GetGene( sLoc.m_iNode, i ).GetName( ) );

	return ( vecstrGenes.size( ) - iSize ); }

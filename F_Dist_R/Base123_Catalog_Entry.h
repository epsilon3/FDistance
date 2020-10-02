// Base123_Catalog_Entry.h

////////////////////////////////////////////////////////////////////////////////
//
//  Base123_Catalog_Entry class (header) encapsulates all the functionality of the Base123 accession catalog
//      entry as drawn from the NCBI .gbk or .dat (influenza, only) file format(s);
//
//  Developed by Stephen Donald Huff, PhD (Stephen.Huff.3@us.af.mil)
//  Biological Informatics Group, RHDJ, 711HPW, United States Air Force Research Laboratory
//  9 July 2016
//  (All Rights Reserved)
//
////////////////////////////////////////////////////////////////////////////////

#pragma once

using namespace std;

class CBase123_Catalog_Entry
{
	//  Initialization

public:

	//  Constructur
	CBase123_Catalog_Entry();
	//  Destructor
	virtual ~CBase123_Catalog_Entry();

	//  Interface (public)

public:

	vector<structCDS>* GetCDSCollectionPointer();
	bool GetCDSCollection(vector<structCDS>& vCDSs);
	bool GetIsSet();
	void SetIsSet(bool bIsSet);
	string GetAccession();
	void SetAccession(string strAccession);
	string GetNameID();
	void SetNameID(string strNameID);
	long GetLength();
	void SetLength(long lLength);
	string GetMoleculeType();
	void SetMoleculeType(string strMoleculeType);
	string GetStrandedness();
	void SetStrandedness(string strStrandedness);
	string GetStrandednessType();
	void SetStrandednessType(string strStrandednessType);
	string GetStrandednessDirection();
	void SetStrandednessDirection(string strStrandednessDirection);
	string GetDNAIntermediate();
	void SetDNAIntermediate(string strDNAIntermediate);
	string GetRNAIntermediate();
	void SetRNAIntermediate(string strDNAIntermediate);
	string GetRetroTranscriptase();
	void SetRetroTranscriptase(string strRetroTranscriptase);
	string GetCompleteness();
	void SetCompleteness(string strCompleteness);
	string GetDate();
	void SetDate(string strDate);
	string GetHost();
	void SetHost(string strHost);
	string GetHostAge();
	void SetHostAge(string strHostAge);
	string GetHostGender();
	void SetHostGender(string strHostGender);
	string GetChromosomeSegment();
	void SetChromosomeSegment(string strChromosomeSegment);
	string GetSeroType();
	void SetSeroType(string strSeroType);
	string GetLocale();
	void SetLocale(string strLocale);
	string GetSatelliteStatus();
	void SetSatelliteStatus(string strSatelliteStatus);
	string GetTranscriptStatus();
	void SetTranscriptStatus(string strTranscriptStatus);
	string GetViralGroup();
	void SetViralGroup(string strViralGroup);
	void CDSAdd(structCDS structStartStop);
	void CDSAdd(string strNameID, long lStart, long lStop, string strCompleteness);
	void CDSRemove(int nIndex);
	void CDSGet(int nIndex, structCDS& stGet);
	structCDS* CDSGet(int nIndex);
	int CDSGetCount();
	string CDSGetIsComplement(int nIndex);
	void CDSSetIsComplement(int nIndex, string strIsComplement);
	string CDSGetNameID(int nIndex);
	void CDSSetNameID(int nIndex, string strNameID);
	long CDSGetStart(int nIndex);
	void CDSSetStart(int nIndex, long lStart);
	long CDSGetStop(int nIndex);
	void CDSSetStop(int nIndex, long lStop);
	string CDSGetCompleteness(int nIndex);
	void CDSSetCompleteness(int nIndex, string strCompleteness);
	void ClearCDS();
	string GetEntry();
	string GetDemographics();
	bool SetEntry(string strEntry);

	//  Implementation (private)
	string GetEntryLine();
	string GetCDSLine();

private:

	//  GBK fiile header format (GenBank flat file header)
	//  NC_010314               1090 bp ss - DNA     circular VRL 20 - OCT - 2015

	//  Dat file header format
	//  LN612606	Human	4	H1N1	Malaysia	2009 / 08 / 11	1696	Influenza A virus(A / Malaysia / 2097724 / 2009(H1N1))			nc

	//  Accession (locus name)
	string m_strAccession;

	//  Name/ID
	string m_strNameID;

	//  Sequence length (bp)
	long m_lLength;

	//  Strandedness
	//  Single [ss]
	//  Double [ds]
	string m_strStrandedness;

	//  Strandedness type
	//  Circular [c]
	//  Linear [l]
	string m_strStrandednessType;

	//  Molecule type
	//  [dna]
	//  [rna]
	string m_strMoleculeType;

	//  Sequence direction
	//  Positive [+]
	//  Negative [-]
	string m_strStrandednessDirection;

	//  DNA intermediate (viruses)
	//  Yes [y]
	//  No or empty [n]
	string m_strDNAIntermediate;

	//  RNA intermediate (viruses)
	//  Yes [y]
	//  No or empty [n]
	string m_strRNAIntermediate;

	//  Retro-transcriptase (viruses)
	//  Yes [y]
	//  No or empty [n]
	string m_strRetroTranscriptase;

	//  The influenza_na.dat and influenza_aa.dat files have an additional field in the last column to indicate the completeness of a sequence
	//  "c" for complete sequences that include start and stop codons;
	//  "nc" for nearly complete sequences that are missing only start and / or stop codons;
	//  "p" for partial sequences.
	//  Complete [c]
	//  Partial [p]
	//  Near [nc]
	string m_strCompleteness;

	//  Date
	string m_strDate;

	//  Host
	string m_strHost;

	//  Host age
	string m_strHostAge;

	//  Host gender
	//  Male [m]
	//  Female [f]
	//  Unknown [u]
	string m_strHostGender;

	//  Segment
	string m_strChromosomeSegment;

	//  Serotype
	string m_strSeroType;

	//  Locale
	string m_strLocale;

	//  Is satellite [y/n]
	string m_strSatelliteStatus;

	//  Is transcript [y/n]
	string m_strTranscriptStatus;

	//  Viral group
	// i, ii , iii, iv, v, vi, vii
	string m_strViralGroup;
	
	//  CDSs for this entry
	vector<structCDS> m_vCDSs;

	//  Entry is set, if true
	bool m_bIsSet;
};


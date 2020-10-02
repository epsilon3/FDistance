// Base123_Catalog.h

////////////////////////////////////////////////////////////////////////////////
//
//  Base123_Catalog class (header) encapsulates all the functionality of the Base123 accession catalog
//      as drawn from the NCBI .gbk or .dat (influenza, only) file format(s);
//
//  Developed by Stephen Donald Huff, PhD (Stephen.Huff.3@us.af.mil)
//  Biological Informatics Group, RHDJ, 711HPW, United States Air Force Research Laboratory
//  9 July 2016
//  (All Rights Reserved)
//
////////////////////////////////////////////////////////////////////////////////

#pragma once

using namespace std;

#include <deque>

class CBase123_Catalog
{
	//  Initialization

public:

	//  Constructor
	CBase123_Catalog(long lMaxSize);
	//  Destructor
	virtual ~CBase123_Catalog();

	//  Interface (public)

public:

	vector<structCDS>* GetCDSCollectionPointerByAccession(string strAccession);
	bool GetCDSCollectionByAccession(string strAccession, vector<structCDS>& vCDSs);
	vector<structCDS>* GetCDSCollectionPointerByIndex(long lIndex);
	bool GetCDSCollectionByIndex(long lIndex, vector<structCDS>& vCDSs);
	bool AddEntry(CBase123_Catalog_Entry& eAdd);
	bool SetEntryAtIndex(CBase123_Catalog_Entry& eSet, long lIndex);
	bool GetEntryByIndex(int nIndex, CBase123_Catalog_Entry& ceGet);
	bool GetEntryByAccession(string strAccession, CBase123_Catalog_Entry& ceGet);
	string GetDemographicsHeader();
	bool RemoveEntry(int nIndex);
	bool OpenCatalog(string strInputFilePathName);
	bool CloseCatalog();
	bool WriteCatalog(string strOutputFilePathName);
	bool CreateDatCatalog(string strDatInputFilePathName, string strNADatInputFilePathName, string strOutputCatalogFilePathName, string strErrorFilePathName, int nMaxProcs);
	bool CreateGBKCatalog(string strInputFilePathNameList, string strInputFilePathNameTransform, string strOutputCatalogFilePathName, string strErrorFilePathName, int nMaxProcs);

	//  Implementation (private)

private:

	bool SetGBKEntryFeatures(CBase123_Catalog_Entry& eSet, string& strFeatures);
	bool CatalogGBKFile(string& strFileText, long lIndex);
	bool CatalogDatFile(string& strDatFileText, string& strNADatFileText, string strErrorFilePathName, int nMaxProcs);
	bool ParseStartsStops(string& strParts, vector<long>& vStartsStops);
	string GetEntryHeader();
	string GetCDSHeader();
	bool ClearEntries();

	//  GBK Features
	string m_strGBKFeature_Locus;
	string m_strGBKFeature_Definition;
	string m_strGBKFeature_Version;
	string m_strGBKFeature_Organism;
	string m_strGBKFeature_Host;
	string m_strGBKFeature_CollectionDate;
	string m_strGBKFeature_Chromosome;
	string m_strGBKFeature_Segment;
	string m_strGBKFeature_Country;
	string m_strGBKFeature_CDS;
	string m_strGBKFeature_Origin;
	string m_strGBKFeature_Gene;
	string m_strGBKFeature_Product;
	string m_strGBKFeature_NextFeature;
	string m_strGBKFeature_GeneID;
	string m_strGBKFeature_Translation;

	//  Catalog entries
	vector<CBase123_Catalog_Entry> m_vEntries;
};


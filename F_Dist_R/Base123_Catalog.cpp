// Base123_Catalog.cpp

////////////////////////////////////////////////////////////////////////////////
//
//  Base123_Catalog class encapsulates all the functionality of the Base123 accession catalog
//      as drawn from the NCBI .gbk or .dat (influenza, only) file format(s);
//
//  Developed by Stephen Donald Huff, PhD (Stephen.Huff.3@us.af.mil)
//  Biological Informatics Group, RHDJ, 711HPW, United States Air Force Research Laboratory
//  9 July 2016
//  (All Rights Reserved)
//
////////////////////////////////////////////////////////////////////////////////

#include "F_Dist_R.h"
#include "Base123_Catalog_Entry.h"
#include "Base123_Catalog.h"
#include "Base123_Utilities.h"

#include <sstream>
#include <omp.h>

//  Initialization

////////////////////////////////////////////////////////////////////////////////
//
//  Constructs the CBase123_Catalog class object
//
////////////////////////////////////////////////////////////////////////////////
//
//  [long] lMaxSize:  maximum size expected for this catalog (indexed upper limit)
//         
////////////////////////////////////////////////////////////////////////////////

CBase123_Catalog::CBase123_Catalog(long lMaxSize)
{
	try
	{
		//  GBK Features MUST be in lowercase for comparisons
		m_strGBKFeature_Locus = "locus       ";
		m_strGBKFeature_Organism = "  organism  ";
		m_strGBKFeature_Definition = "definition  ";
		m_strGBKFeature_Version = "version     ";
		m_strGBKFeature_Host = "                     /host=\"";
		m_strGBKFeature_Chromosome = "                     /chromosome=\"";
		m_strGBKFeature_Segment = "                     /segment=\"";
		m_strGBKFeature_Country = "                     /country=\"";
		m_strGBKFeature_CollectionDate = "                     /collection_date=\"";
		m_strGBKFeature_CDS = "     cds             ";
		m_strGBKFeature_Gene = "                     /gene=\"";
		m_strGBKFeature_Product = "                     /product=\"";
		m_strGBKFeature_NextFeature = "                     /";
		m_strGBKFeature_GeneID = "                     /db_xref=\"geneid:";
		m_strGBKFeature_Translation = "                     /translation=\"";
		m_strGBKFeature_Origin = "origin";

		m_vEntries.clear();
		m_vEntries.resize(lMaxSize);
	}
	catch (exception ex)
	{
		cout << "ERROR [CBase123_Catalog] Exception Code:  " << ex.what() << "\n";
	}
}

////////////////////////////////////////////////////////////////////////////////
//
//  Destructs the CBase123_Catalog class object
//
////////////////////////////////////////////////////////////////////////////////
//
//
//         
////////////////////////////////////////////////////////////////////////////////

CBase123_Catalog::~CBase123_Catalog()
{
	try
	{
	}
	catch (exception ex)
	{
		cout << "ERROR [~CBase123_Catalog] Exception Code:  " << ex.what() << "\n";
	}
}

//  Interface (public)

////////////////////////////////////////////////////////////////////////////////
//
//  Gets a pointer to the indicated CDS colletion vector<structCDS> by accession search
//
////////////////////////////////////////////////////////////////////////////////
//
//  [string] strAccession:  accession to search and return
//                      :  returns  a pointer to the requested CDS collection vector<structCDS>, if successful; else, NULL
//         
////////////////////////////////////////////////////////////////////////////////

vector<structCDS>* CBase123_Catalog::GetCDSCollectionPointerByAccession(string strAccession)
{
	try
	{
		if (!strAccession.empty())
		{
			for (long lCount = 0; lCount < m_vEntries.size(); lCount++)
			{
				if (m_vEntries[lCount].GetAccession() == strAccession)
					return m_vEntries[lCount].GetCDSCollectionPointer();
			}
		}
		else
		{
			ReportTimeStamp("[GetCDSCollectionPointerByAccession]", "ERROR:  Accession is Empty");
		}
	}
	catch (exception ex)
	{
		cout << "ERROR [GetCDSCollectionPointerByAccession] Exception Code:  " << ex.what() << "\n";
	}

	return NULL;
}

////////////////////////////////////////////////////////////////////////////////
//
//  Gets the indicated CDS colletion vector<structCDS> by accession search
//
////////////////////////////////////////////////////////////////////////////////
//
//  [string] strAccession :  accession to search and return
//  [vector<structCDS>&] vCDSs:  accession to search and return
//                       :  returns true, if successful; else, false
//         
////////////////////////////////////////////////////////////////////////////////

bool CBase123_Catalog::GetCDSCollectionByAccession(string strAccession, vector<structCDS>& vCDSs)
{
	try
	{
		if (!strAccession.empty())
		{
			for (long lCount = 0; lCount < m_vEntries.size(); lCount++)
			{
				if (m_vEntries[lCount].GetAccession() == strAccession)
					return m_vEntries[lCount].GetCDSCollection(vCDSs);
			}
		}
		else
		{
			ReportTimeStamp("[GetCDSCollectionByAccession]", "ERROR:  Accession is Empty");
		}
	}
	catch (exception ex)
	{
		cout << "ERROR [GetCDSCollectionByAccession] Exception Code:  " << ex.what() << "\n";
	}

	return false;
}

////////////////////////////////////////////////////////////////////////////////
//
//  Gets a pointer to the indicated CDS colletion vector<structCDS> by accession search
//
////////////////////////////////////////////////////////////////////////////////
//
//  [long] lINdex:  index to return
//              :  returns  a pointer to the requested CDS collection vector<structCDS>, if successful; else, false
//         
////////////////////////////////////////////////////////////////////////////////

vector<structCDS>* CBase123_Catalog::GetCDSCollectionPointerByIndex(long lIndex)
{
	try
	{
		if ((lIndex >= 0) && (lIndex < m_vEntries.size()))
		{
			return m_vEntries[lIndex].GetCDSCollectionPointer();
		}
		else
		{
			ReportTimeStamp("[GetCDSCollectionPointerByIndex]", "ERROR:  Index [" + ConvertLongToString(lIndex) + "] is Out of Range [0:" + ConvertLongToString((long)m_vEntries.size() - 1) + "]");
		}
	}
	catch (exception ex)
	{
		cout << "ERROR [GetCDSCollectionPointerByIndex] Exception Code:  " << ex.what() << "\n";
	}

	return NULL;
}

////////////////////////////////////////////////////////////////////////////////
//
//  Gets the indicated CDS colletion vector<structCDS> by index search
//
////////////////////////////////////////////////////////////////////////////////
//
//  [long] lIndex         :  index to return
//  [vector<structCDS>&] vCDSs:  accession to search and return
//                       :  returns true, if successful; else, false
//         
////////////////////////////////////////////////////////////////////////////////

bool CBase123_Catalog::GetCDSCollectionByIndex(long lIndex, vector<structCDS>& vCDSs)
{
	try
	{
		if ((lIndex >= 0) && (lIndex < m_vEntries.size()))
		{
			return m_vEntries[lIndex].GetCDSCollection(vCDSs);
		}
		else
		{
			ReportTimeStamp("[GetCDSCollectionByIndex]", "ERROR:  Index [" + ConvertLongToString(lIndex) + "] is Out of Range [0:" + ConvertLongToString((long)m_vEntries.size() - 1) + "]");
		}
	}
	catch (exception ex)
	{
		cout << "ERROR [GetCDSCollectionByIndex] Exception Code:  " << ex.what() << "\n";
	}

	return false;
}

////////////////////////////////////////////////////////////////////////////////
//
//  Sets a catalog entry at a given index
//
////////////////////////////////////////////////////////////////////////////////
//
//  [CBase123_Catalog_Entry] eSet:  catalog entry to set
//  [long] lIndex                :  index position to set
//                              :  returns true, if successful; else, false
//         
////////////////////////////////////////////////////////////////////////////////

bool CBase123_Catalog::SetEntryAtIndex(CBase123_Catalog_Entry& eSet, long lIndex)
{
	try
	{
		if ((lIndex >= 0) && (lIndex < m_vEntries.size()))
		{
			eSet.SetIsSet(true);

			m_vEntries[lIndex] = eSet;

			return true;
		}
		else
		{
			ReportTimeStamp("[SetEntryAtIndex]", "ERROR:  Index [" + ConvertLongToString(lIndex) + "] is Out of Range [0:" + ConvertLongToString((long)m_vEntries.size() - 1) + "]");
		}
	}
	catch (exception ex)
	{
		cout << "ERROR [SetEntryAtIndex] Exception Code:  " << ex.what() << "\n";
	}

	return false;
}

////////////////////////////////////////////////////////////////////////////////
//
//  Adds a catalog entry
//
////////////////////////////////////////////////////////////////////////////////
//
//  [CBase123_Catalog_Entry] eAdd:  catalog entry to add
//                              :  returns true, if successful; else, false
//         
////////////////////////////////////////////////////////////////////////////////

bool CBase123_Catalog::AddEntry(CBase123_Catalog_Entry& eAdd)
{
	try
	{
		eAdd.SetIsSet(true);

		m_vEntries.push_back(eAdd);

		return true;
	}
	catch (exception ex)
	{
		cout << "ERROR [AddEntry] Exception Code:  " << ex.what() << "\n";
	}

	return false;
}

////////////////////////////////////////////////////////////////////////////////
//
//  Gets the catalog entry at the indicated index
//
////////////////////////////////////////////////////////////////////////////////
//
//  [int] nIndex                   :  index to get
//  [CBase123_Catalog_Entry&] ceGet:  catalog entry to return
//                                :  returns true, if successful; else, false
//         
////////////////////////////////////////////////////////////////////////////////

bool CBase123_Catalog::GetEntryByIndex(int nIndex, CBase123_Catalog_Entry& ceGet)
{
	try
	{
		//  If the index is in range
		if ((nIndex >= 0) && (nIndex < m_vEntries.size()))
		{
			ceGet = m_vEntries[nIndex];

			return true;
		}
		else
		{
			ReportTimeStamp("[GetEntryByIndex]", "ERROR:  Catalog Index [" + ConvertIntToString(nIndex) + "] is Out of Range [0:" + ConvertIntToString((int)m_vEntries.size() - 1) + "]");
		}
	}
	catch (exception ex)
	{
		cout << "ERROR [GetEntryByIndex] Exception Code:  " << ex.what() << "\n";
	}

	return false;
}

////////////////////////////////////////////////////////////////////////////////
//
//  Gets the catalog entry for indicated accession
//
////////////////////////////////////////////////////////////////////////////////
//
//  [string] strAccession          :  accession to get
//  [CBase123_Catalog_Entry&] ceGet:  catalog entry to return
//                                :  returns true, if successful; else, false
//         
////////////////////////////////////////////////////////////////////////////////

bool CBase123_Catalog::GetEntryByAccession(string strAccession, CBase123_Catalog_Entry& ceGet)
{
	try
	{
		//  If the accession is not empty
		if (!strAccession.empty())
		{
			for (int nCount = 0; nCount < m_vEntries.size(); nCount++)
			{
				if (m_vEntries[nCount].GetAccession() == strAccession)
				{
					ceGet = m_vEntries[nCount];

					return true;
				}

				if (!m_vEntries[nCount].GetIsSet())
					break;
			}
		}
		else
		{
			ReportTimeStamp("[GetEntryByAccession]", "ERROR:  Accession is Empty");
		}
	}
	catch (exception ex)
	{
		cout << "ERROR [GetEntryByAccession] Exception Code:  " << ex.what() << "\n";
	}

	return false;
}

////////////////////////////////////////////////////////////////////////////////
//
//  Gets the demographics (entry) header, only
//
////////////////////////////////////////////////////////////////////////////////
//
//  :  returns demographics header, if successful; else, an empty string
//         
////////////////////////////////////////////////////////////////////////////////

string CBase123_Catalog::GetDemographicsHeader()
{
	try
	{
		return GetEntryHeader();
	}
	catch (exception ex)
	{
		cout << "ERROR [GetDemographicsHeader] Exception Code:  " << ex.what() << "\n";
	}

	return "";
}

////////////////////////////////////////////////////////////////////////////////
//
//  Removes the catalog entry at the indicated index
//
////////////////////////////////////////////////////////////////////////////////
//
//  [int] nIndex:  index to remove
//             :  returns true, if successful; else, false
//         
////////////////////////////////////////////////////////////////////////////////

bool CBase123_Catalog::RemoveEntry(int nIndex)
{
	try
	{
		//  If the index is in range
		if ((nIndex >= 0) && (nIndex < m_vEntries.size()))
		{
			m_vEntries.erase(m_vEntries.begin() + nIndex);

			return true;
		}
		else
		{
			ReportTimeStamp("[RemoveEntry]", "ERROR:  Catalog Index [" + ConvertIntToString(nIndex) + "] is Out of Range [0:" + ConvertIntToString((int)m_vEntries.size() - 1) + "]");
		}
	}
	catch (exception ex)
	{
		cout << "ERROR [RemoveEntry] Exception Code:  " << ex.what() << "\n";
	}

	return false;
}

////////////////////////////////////////////////////////////////////////////////
//
//  Opens and creates the catalog from file
//
////////////////////////////////////////////////////////////////////////////////
//
//  [string] strInputFilePathName:  file path name of the catalog
//                              :  returns true, if successful; else, false
//         
////////////////////////////////////////////////////////////////////////////////

bool CBase123_Catalog::OpenCatalog(string strInputFilePathName)
{
	//  The file text to open
	string strInputFileText = "";
	//  The file line vector
	vector<string> vLines;

	try
	{
		//  If the file path name is not empty
		if (!strInputFilePathName.empty())
		{
			//  Get the file text
			GetFileText(strInputFilePathName, strInputFileText);

			//  If file text is not empty
			if (!strInputFileText.empty())
			{
				//  SPlit the file into lines
				SplitString(strInputFileText, '\n', vLines);

				//  If lines are filled
				if (vLines.size() > 0)
				{
					//  Iterate through the lines, create catalog entries; 1, to skip header
					for (long lCount = 1; lCount < vLines.size(); lCount++)
					{
						//  If line is not empty
						if (!vLines[lCount].empty())
						{
							CBase123_Catalog_Entry eAdd;

							eAdd.SetEntry(vLines[lCount]);

							if ((lCount - 1 >= 0) && (lCount - 1 < m_vEntries.size()))
								SetEntryAtIndex(eAdd, lCount - 1);
							else
							{
								ReportTimeStamp("[OpenCatalog]", "ERROR:  Catalog Entry Index [" + ConvertLongToString(lCount) + "] Exceeds Index Range [0:" + ConvertLongToString((long)m_vEntries.size()) + "]");

								return false;
							}
						}
						//  No entry, there should not be empty lines
					}

					return true;
				}
				else
				{
					ReportTimeStamp("[OpenCatalog]", "ERROR:  Catalog File [" + strInputFilePathName + "] Container is Not Set");
				}
			}
			else
			{
				ReportTimeStamp("[OpenCatalog]", "ERROR:  Catalog File [" + strInputFilePathName + "] Text is Empty");
			}
		}
		else
		{
			ReportTimeStamp("[OpenCatalog]", "ERROR:  Catalog File Path Name is Empty");
		}
	}
	catch (exception ex)
	{
		cout << "ERROR [OpenCatalog] Exception Code:  " << ex.what() << "\n";
	}

	return false;
}

////////////////////////////////////////////////////////////////////////////////
//
//  Closes the catalog and clears datat
//
////////////////////////////////////////////////////////////////////////////////
//
//  :  returns true, if successful; else, false
//         
////////////////////////////////////////////////////////////////////////////////

bool CBase123_Catalog::CloseCatalog()
{
	try
	{
		//  Catalog entries
		ClearEntries();
	}
	catch (exception ex)
	{
		cout << "ERROR [CloseCatalog] Exception Code:  " << ex.what() << "\n";
	}

	return false;
}

////////////////////////////////////////////////////////////////////////////////
//
//  Writes the catalog to file
//
////////////////////////////////////////////////////////////////////////////////
//
//  [string] strOutputFilePathName:  file path name of the catalog
//                               :  returns true, if successful; else, false
//         
////////////////////////////////////////////////////////////////////////////////

bool CBase123_Catalog::WriteCatalog(string strOutputFilePathName)
{
	//  The file text to write
	string strOutputFileText = "";

	try
	{
		//  If the file path name is not empty
		if (!strOutputFilePathName.empty())
		{
			//  Write headers
			strOutputFileText = GetEntryHeader() + "|" + GetCDSHeader() + "\n" + strOutputFileText;

			//  Iterate entries to fill file text
			for (long lCount = 0; lCount < m_vEntries.size() ; lCount++)
			{
				//  If entry is set
				if (m_vEntries[lCount].GetIsSet())
				{
					//  Append new line to second+ lines
					if (lCount > 0)
						strOutputFileText += '\n';

					//  Concatenate the file text
					strOutputFileText += m_vEntries[lCount].GetEntry();
				}
			}

			//  Write the text to file
			if (!WriteFileText(strOutputFilePathName, strOutputFileText))
				ReportTimeStamp("[WriteCatalog]", "ERROR:  Catalog File [" + strOutputFilePathName + "] Write Failed");
			else
				return true;
		}
		else
		{
			ReportTimeStamp("[WriteCatalog]", "ERROR:  Catalog File Path Name is Empty");
		}
	}
	catch (exception ex)
	{
		cout << "ERROR [WriteCatalog] Exception Code:  " << ex.what() << "\n";
	}

	return false;
}

////////////////////////////////////////////////////////////////////////////////
//
//  Creates a BIG format genome data catalog from NCBI .dat/_na.dat format
//         
////////////////////////////////////////////////////////////////////////////////
//
//  [string] strDatInputFilePathName     :  input .dat file path name
//  [string] strNADatInputFilePathName   :  input _na.dat file path name
//  [string] strOutputCatalogFilePathName:  output catalog file path name  
//  [string] strErrorFilePathName        :  error file path name
//  [int] nMaxProcs                      :  maximum processors for openMP
//                                      :  returns true, if successful; else, false  
//         
////////////////////////////////////////////////////////////////////////////////

bool CBase123_Catalog::CreateDatCatalog(string strDatInputFilePathName, string strNADatInputFilePathName, string strOutputCatalogFilePathName, string strErrorFilePathName, int nMaxProcs)
{
	//  .dat file text
	string strDatFileText = "";
	//  _na.dat file text
	string strNADatFileText = "";

	try
	{
		//  If .dat file path name is not empty
		if (!strDatInputFilePathName.empty())
		{
			//  If _na.dat file path name is not empty
			if (!strNADatInputFilePathName.empty())
			{
				//  If DAT catalog file path name is not empty
				if (!strOutputCatalogFilePathName.empty())
				{
					//  Read .dat file text
					if (GetFileText(strDatInputFilePathName, strDatFileText))
					{
						//  Read _na.dat file text
						if (GetFileText(strNADatInputFilePathName, strNADatFileText))
						{
							//  Catalog the data
							if (CatalogDatFile(strDatFileText, strNADatFileText, strErrorFilePathName, nMaxProcs))
							{
								//  Write the catalog
								if (!WriteCatalog(strOutputCatalogFilePathName))
									ReportTimeStamp("[CreateDatCatalog]", "ERROR:  _na.dat Catalog [" + strOutputCatalogFilePathName + "] Write Failed");
								else
									return true;
							}
							else
							{
								ReportTimeStamp("[CreateDatCatalog]", "ERROR:  _na.dat Catalog [" + strOutputCatalogFilePathName + "] / [" + strDatInputFilePathName + "] / [" + strNADatInputFilePathName + "] Creation Failed");
							}
						}
						else
						{
							ReportTimeStamp("[CreateDatCatalog]", "ERROR:  _na.dat Input File [" + strNADatInputFilePathName + "] Open Failed");
						}
					}
					else
					{
						ReportTimeStamp("[CreateDatCatalog]", "ERROR:  .dat Input File [" + strDatInputFilePathName + "] Open Failed");
					}
				}
				else
				{
					ReportTimeStamp("[CreateDatCatalog]", "ERROR:  DAT Catalog File Path Name is Empty");
				}
			}
			else
			{
				ReportTimeStamp("[CreateDatCatalog]", "ERROR:  _na.dat Input File Path Name is Empty");
			}
		}
		else
		{
			ReportTimeStamp("[CreateDatCatalog]", "ERROR:  .dat Input File Path Name is Empty");
		}
	}
	catch (exception ex)
	{
		cout << "ERROR [CreateDatCatalog] Exception Code:  " << ex.what() << "\n";
	}

	return false;
}

////////////////////////////////////////////////////////////////////////////////
//
//  Creates a BIG format genome data catalog from NCBI .gbk format
//         
////////////////////////////////////////////////////////////////////////////////
//
//  [string] strInputFilePathNameList     :  input file path name list (in NCBI .gbk format)
//  [string] strInputFilePathNameTransform:  input file path name transform (includes string replacements, see help)
//  [string] strOutputCatalogFilePathName :  output catalog file path name  
//  [string] strErrorFilePathName         :  error file path name
//  [int] nMaxProcs                       :  maximum processors for openMP
//                                       :  returns true, if successful; else, false  
//         
////////////////////////////////////////////////////////////////////////////////

bool CBase123_Catalog::CreateGBKCatalog(string strInputFilePathNameList, string strInputFilePathNameTransform, string strOutputCatalogFilePathName, string strErrorFilePathName, int nMaxProcs)
{
	//  List lock
	omp_lock_t lockList;
	//  Input list file text
	string strInputListFileText = "";
	//  File path name vector<string>
	vector<string> vFilePathNames;
	//  Error list file entry text
	vector<string> vErrorEntries;
	//  Error file text
	string strErrorFileText = "";

	try
	{
		//  If input file path name
		if (!strInputFilePathNameList.empty())
		{
			//  If output catalog file path name is not empty
			if (!strOutputCatalogFilePathName.empty())
			{
				//  Get .gbk list file text
				if (GetFileText(strInputFilePathNameList, strInputListFileText))
				{
					//  Split list
					SplitString(strInputListFileText, '\n', vFilePathNames);

					//  If list contains entries
					if (vFilePathNames.size() > 0)
					{
						//  Initialize errof file vector<string>
						vErrorEntries.resize(vFilePathNames.size());

						//  Initialize time stamp lock
						omp_init_lock(&lockList);

						//  Declare omp parallel
						#pragma omp parallel num_threads(nMaxProcs)
						{
							//  omp loop
							#pragma omp for
							for (long lCount = 0; lCount < vFilePathNames.size(); lCount++)
							{
								//  Test max procs
								if (lCount == 0)
								{
									omp_set_lock(&lockList);
									ReportTimeStamp("[CatalogGBKFile]", "NOTE:  Thread Count = " + ConvertIntToString(omp_get_num_threads()));
									omp_unset_lock(&lockList);
								}

								if (lCount % 10000 == 0)
								{
									omp_set_lock(&lockList);
									ReportTimeStamp("[CatalogGBKFile]", "NOTE:  Processing Entry [" + ConvertLongToString(lCount) + "] [" + vFilePathNames[lCount] + "]");
									omp_unset_lock(&lockList);
								}

								//  If file path name is not empty
								if (!vFilePathNames[lCount].empty())
								{
									//  GBK file text
									string strGBKFileText = "";
									//  Working file path name
									string strWorkingFilePathName = "";

									//  If input file path name transform is not empty
									if (!strInputFilePathNameTransform.empty())
										strWorkingFilePathName = TransformFilePathName(vFilePathNames[lCount], strInputFilePathNameTransform, "");
									else
										strWorkingFilePathName = vFilePathNames[lCount];

									//  Read .gbk file text
									if (GetFileText(strWorkingFilePathName, strGBKFileText))
									{
										//  Catalog the data
										if (!CatalogGBKFile(strGBKFileText, lCount))
										{
											vErrorEntries[lCount] = strWorkingFilePathName + "~Catalog Failed\n";

											omp_set_lock(&lockList);
											ReportTimeStamp("[CreateGBKCatalog]", "ERROR:  .gbk Catalog Entry [" + strWorkingFilePathName + "] Failed");
											omp_unset_lock(&lockList);
										}
									}
									else
									{
										vErrorEntries[lCount] = strWorkingFilePathName + "~Open Failed\n";

										omp_set_lock(&lockList);
										ReportTimeStamp("[CreateGBKCatalog]", "ERROR:  .gbk Input File [" + strWorkingFilePathName + "] Open Failed");
										omp_unset_lock(&lockList);
									}
								}
							}
						}

						//  Destroy time stamp lock
						omp_destroy_lock(&lockList);

						//  Write error file
						if (!strErrorFilePathName.empty())
						{
							//  Iterate list file vector and concatenate file
							for (long lCount = 0; lCount < vErrorEntries.size(); lCount++)
							{
								if (!vErrorEntries[lCount].empty())
									strErrorFileText += vErrorEntries[lCount];
							}

							//  Write the file
							if (!WriteFileText(strErrorFilePathName, strErrorFileText))
								ReportTimeStamp("[CreateGBKCatalog]", "ERROR:  GBK Catalog Error File [" + strErrorFilePathName + "] Write Failed");
						}

						//  Write the catalog
						if (WriteCatalog(strOutputCatalogFilePathName))
							return true;
						else
						{
							ReportTimeStamp("[CreateGBKCatalog]", "ERROR:  .gbk Catalog File [" + strOutputCatalogFilePathName + "] Write Failed");
						}
					}
					else
					{
						ReportTimeStamp("[CreateGBKCatalog]", "ERROR:  .gbk File Path Name List Container is Not Set");
					}
				}
				else
				{
					ReportTimeStamp("[CreateGBKCatalog]", "ERROR:  .gbk List File Text Open Failed");
				}
			}
			else
			{
				ReportTimeStamp("[CreateGBKCatalog]", "ERROR:  .gbk Catalog File Path Name is Empty");
			}
			}
		else
		{
			ReportTimeStamp("[CreateGBKCatalog]", "ERROR:  .gbk Input File Path Name List is Empty");
		}
	}
	catch (exception ex)
	{
		cout << "ERROR [CreateGBKCatalog] Exception Code:  " << ex.what() << "\n";
	}

	return false;
}

//  Implementation (private)

////////////////////////////////////////////////////////////////////////////////
//
//  Sets a catalog entry features based on a .gbk catalog feature line
//
////////////////////////////////////////////////////////////////////////////////
//
//  [CBase123_Catalog_Entry&] eSet:  catalog entry to set
//  [string&] strFeatures         :  .gbk catalog feature line to parse (must be lowercase)
//                               :  returns true, if successful; else, false
//         
////////////////////////////////////////////////////////////////////////////////

bool CBase123_Catalog::SetGBKEntryFeatures(CBase123_Catalog_Entry& eSet, string& strFeatures)
{
	//  Locus fields
	vector<string> vFields;
	try
	{
		if (!strFeatures.empty())
		{
			//  Split the locus; all double spaces create empty entries, which are not pushed onto the vector
			SplitString(strFeatures, ' ', vFields);

			//  Parse and set features, where present
			for (int nCountFields = 0; nCountFields < vFields.size(); nCountFields++)
			{
				//  Set strandedness and molecule type
				if ((vFields[nCountFields] == "dsdna") || (vFields[nCountFields] == "ds-dna"))
				{
					eSet.SetStrandedness("ds");
					eSet.SetMoleculeType("dna");
					eSet.SetStrandednessDirection("+/-");
				}
				else if ((vFields[nCountFields] == "ssdna") || (vFields[nCountFields] == "ss-dna"))
				{
					eSet.SetStrandedness("ss");
					eSet.SetMoleculeType("dna");
					//  Explicitly clear direction, in case the solo [dna/rna] case set it previously
					eSet.SetStrandednessDirection("");
				}
				else if ((vFields[nCountFields] == "dsrna") || (vFields[nCountFields] == "ds-rna"))
				{
					eSet.SetStrandedness("ds");
					eSet.SetMoleculeType("rna");
					eSet.SetStrandednessDirection("+/-");
				}
				else if ((vFields[nCountFields] == "ssrna") || (vFields[nCountFields] == "ss-rna"))
				{
					eSet.SetStrandedness("ss");
					eSet.SetMoleculeType("rna");
					//  Explicitly clear direction, in case the solo [dna/rna] case set it previously
					eSet.SetStrandednessDirection("");
				}
				else if ((vFields[nCountFields] == "mrna") || (vFields[nCountFields] == "m-rna"))
				{
					eSet.SetMoleculeType("rna");
					eSet.SetStrandedness("ss");
					eSet.SetStrandednessDirection("-");
					eSet.SetStrandednessType("l");
					eSet.SetTranscriptStatus("y");
				}

				//  Set strandedness direction
				if ((vFields[nCountFields] == "positive") || (vFields[nCountFields] == "positive-strand"))
					eSet.SetStrandednessDirection("+");
				else if ((vFields[nCountFields] == "negative") || (vFields[nCountFields] == "negative-strand"))
					eSet.SetStrandednessDirection("-");

				//  Strandedness type
				if (vFields[nCountFields] == "linear")
					eSet.SetStrandednessType("l");
				else if (vFields[nCountFields] == "circular")
					eSet.SetStrandednessType("c");

				//  Completeness
				if (vFields[nCountFields] == "complete")
					eSet.SetCompleteness("c");
				else if (vFields[nCountFields] == "partial")
					eSet.SetCompleteness("p");

				//  Set retro-transcriptase and DNA intermediate
				if ((vFields[nCountFields] == "retro") || (vFields[nCountFields] == "retro-transcribing"))
					eSet.SetRetroTranscriptase("y");
			}

			//  Completeness, nearly complete
			if ((strFeatures.find("double stranded") != string::npos) || (strFeatures.find("double-stranded") != string::npos))
			{
				eSet.SetStrandedness("ds");
				eSet.SetStrandednessDirection("+/-");
			}
			else if ((strFeatures.find("single stranded") != string::npos) || (strFeatures.find("single-stranded") != string::npos))
				eSet.SetStrandedness("ss");

			//  Completeness, nearly complete
			if (strFeatures.find("nearly complete") != string::npos)
				eSet.SetCompleteness("nc");

			//  Set DNA intermediate
			if ((strFeatures.find("no dna stage") != string::npos) || (strFeatures.find("no dna virus") != string::npos))
				eSet.SetDNAIntermediate("n");

			//  Set RNA intermediate
			if (strFeatures.find("no rna stage") != string::npos)
				eSet.SetRNAIntermediate("n");

			//  Set retro-transcriptase
			if (strFeatures.find("retro") != string::npos)
				eSet.SetRetroTranscriptase("y");

			//  Set satellite status
			if (strFeatures.find("satellite") != string::npos)
				eSet.SetSatelliteStatus("y");
		}
	}
	catch (exception ex)
	{
		cout << "ERROR [SetGBKEntryFeatures] Exception Code:  " << ex.what() << "\n";
	}

	return false;
}

////////////////////////////////////////////////////////////////////////////////
//
//  Catalogs an influenza .dat file format
//
////////////////////////////////////////////////////////////////////////////////
//
//  [string&] strDatFileText     :  .dat file text to catalog
//  [string&] strNADatFileText   :  _na.dat file text to catalog
//  [string] strErrorFilePathName:  error file path name
//  [int] nMaxProcs              :  max processors to use with openMP
//                              :  returns true, if successful; else, false
//         
////////////////////////////////////////////////////////////////////////////////

bool CBase123_Catalog::CatalogDatFile(string& strDatFileText, string& strNADatFileText, string strErrorFilePathName, int nMaxProcs)
{
	//  List lock
	omp_lock_t lockList;
	//  Catalog entries, vector<string>
	vector<string> vCatalogEntries;
	//  Entry CDSs, vector<string>
	vector<string> vCDSs;
	//  Error file entries
	vector<string> vErrorEntries;
	//  Error file text
	string strErrorFileText = "";

	try
	{
		//  Set min procs
		if (nMaxProcs <= 0)
			nMaxProcs = 1;

		//  If the .dat file text is not empty
		if (!strDatFileText.empty())
		{
			//  If the _na.dat file text is not empty
			if (!strNADatFileText.empty())
			{
				//  .dat file containes CDSs
				SplitString(strDatFileText, '\n', vCDSs);

				//  If CDSs is set
				if (vCDSs.size() > 0)
				{
					//  _na.dat file containes entries
					SplitString(strNADatFileText, '\n', vCatalogEntries);

					//  If entries is set
					if (vCatalogEntries.size() > 0)
					{
						//  Initialize errof file vector<string>
						vErrorEntries.resize(vCatalogEntries.size());

						//  Initialize time stamp lock
						omp_init_lock(&lockList);

						//  Declare omp parallel
						#pragma omp parallel num_threads(nMaxProcs)
						{
							//  omp loop
							#pragma omp for
							for (long lCountEntries = 0; lCountEntries < vCatalogEntries.size(); lCountEntries++)
							{
								//  Test max procs
								if (lCountEntries == 0)
								{
									omp_set_lock(&lockList);
									ReportTimeStamp("[CatalogDatFile]", "NOTE:  Thread Count = " + ConvertIntToString(omp_get_num_threads()));
									omp_unset_lock(&lockList);
								}

								if (lCountEntries % 10000 == 0)
								{
									omp_set_lock(&lockList);
									ReportTimeStamp("[CatalogDatFile]", "NOTE:  Processing Entry [" + ConvertLongToString(lCountEntries) + "]");
									omp_unset_lock(&lockList);
								}

								//  Entry vector<string>
								vector<string> vEntry;

								//  Split the entry, allow empty entries
								SplitStringAllowEmptyEntries(vCatalogEntries[lCountEntries], "\t", vEntry);

								//CY026046	Human	4	H3N2	Austria	2004/01/09	987	Influenza A virus (A/Austria/135553/2004(H3N2))	12	M	p
								//CY026046.Human.4.H3N2.Austria.2004/01/09.987.Influenza A virus (A/Austria/135553/2004(H3N2)).12.M.p
								//  If entry is properly formatted
								if ((vEntry.size() >= 1) && (!vEntry[0].empty()))
								{
									//  Entry to  add
									CBase123_Catalog_Entry eAdd;

									//  Entry accession
									string strAccession = "";
									//  Entry length
									long lLength = 0;

									eAdd.SetAccession(vEntry[0]);
									strAccession = vEntry[0];

									if (vEntry.size() >= 2)
										eAdd.SetHost(vEntry[1]);

									if (vEntry.size() >= 3)
										eAdd.SetChromosomeSegment(vEntry[2]);

									if (vEntry.size() >= 4)
										eAdd.SetSeroType(vEntry[3]);

									if (vEntry.size() >= 5)
										eAdd.SetLocale(vEntry[4]);

									if (vEntry.size() >= 6)
										eAdd.SetDate(vEntry[5]);

									if (vEntry.size() >= 7)
									{
										stringstream(vEntry[6]) >> lLength;
										eAdd.SetLength(lLength);
									}

									if (vEntry.size() >= 8)
										eAdd.SetNameID(vEntry[7]);

									if (vEntry.size() >= 9)
										eAdd.SetHostAge(vEntry[8]);

									if (vEntry.size() >= 10)
										eAdd.SetHostGender(vEntry[9]);

									if (vEntry.size() >= 11)
										eAdd.SetCompleteness(vEntry[10]);

									eAdd.SetMoleculeType("rna");
									eAdd.SetStrandedness("ss");
									eAdd.SetStrandednessDirection("-");
									eAdd.SetStrandednessType("l");
									eAdd.SetViralGroup("v");

									//  Search for and process CDSs for this accession/entry
									for (int nCountCDSs = 0; nCountCDSs < vCDSs.size(); nCountCDSs++)
									{
										//  Parse CDSs for this accession
										if (vCDSs[nCountCDSs].substr(0, strAccession.length()) == strAccession)
										{
											//AB000728	BAA75881	gb|AB000728:4-1128	BAA75882	(gb|AB000728:4-731, 960)	BAA75883	gb|AB000728:709-1128
											//  CDS parts vector
											vector<string> vCDSParts;

											//  Spit the CDS into parts
											SplitString(vCDSs[nCountCDSs], '\t', vCDSParts);

											//  If there is at least one CDS, first part is the accession
											if (vCDSParts.size() >= 3)
											{
												//  Iterate through CDS part(s), each part should be in a pair
												for (int nCountPart = 1; nCountPart < vCDSParts.size() - 1; nCountPart += 2)
												{
													//  CDS start/stop parts
													vector<string> vCDSStartStop;

													//  CDS to add
													string strCDSNameID = "";
													long lCDSStart = 0;
													long lCDSStop = 0;
													string strCDSCompleteness = "";

													//  Name/ID
													strCDSNameID = vCDSParts[nCountPart];

													//  Split start/stop parts
													SplitString(vCDSParts[nCountPart + 1], ':', vCDSStartStop);

													//  If start/stop part is properly Formatted
													if (vCDSStartStop.size() == 2)
													{
														//  starts and stops, section(s)
														vector<string> vStartsStopsSections;

														//  Completeness
														strCDSCompleteness = vCDSStartStop[1];

														//  Split part into numeric section(s)
														SplitString(vCDSStartStop[1], ',', vStartsStopsSections);

														//  Iterate sections, convert to numerics
														for (int nCountSections = 0; nCountSections < vStartsStopsSections.size(); nCountSections++)
														{
															//  starts and stops, pairwise
															vector<long> vStartsStops;

															//  Parse starts and stops into numeric pairs
															ParseStartsStops(vStartsStopsSections[nCountSections], vStartsStops);

															//  Single value, only
															if (vStartsStops.size() == 1)
															{
																//  Set start/stop
																lCDSStart = vStartsStops[0];
																lCDSStop = vStartsStops[0];

																//  Add CDS to entry
																if ((!strCDSNameID.empty()) && (lCDSStart > 0) && (lCDSStop > 0))
																	eAdd.CDSAdd(strCDSNameID, lCDSStart, lCDSStop, strCDSCompleteness);
															}
															//  Paired value
															if (vStartsStops.size() == 2)
															{
																//  Set start/stop
																lCDSStart = vStartsStops[0];
																lCDSStop = vStartsStops[1];

																//  Add CDS to entry
																if ((!strCDSNameID.empty()) && (lCDSStart > 0) && (lCDSStop > 0))
																	eAdd.CDSAdd(strCDSNameID, lCDSStart, lCDSStop, strCDSCompleteness);
															}
														}
													}
													else
													{
														vErrorEntries[lCountEntries] = "Entry [" + ConvertLongToString(lCountEntries) + "] [" + strAccession + "]~CDS Format Error~" + vCDSParts[nCountPart + 1] + "\n";

														omp_set_lock(&lockList);
														ReportTimeStamp("[CatalogDatFile]", "ERROR:  _na.dat Entry [" + vCatalogEntries[lCountEntries] + "] CDS Start/Stop Part [" + vCDSParts[nCountPart + 1] + "] is Not Properly Formatted");
														omp_unset_lock(&lockList);
													}
												}
											}

											//  End search
											break;
										}
									}

									//  If entry has no CDS defined, the entire length is the CDS
									if (eAdd.CDSGetCount() <= 0)
										eAdd.CDSAdd("Molecular Polymer", 1, eAdd.GetLength(), "1.." + ConvertLongToString(eAdd.GetLength()));

									//  Add the entry to the catalog
									SetEntryAtIndex(eAdd, lCountEntries);
								}
								else
								{
									vErrorEntries[lCountEntries] = "Entry [" + ConvertLongToString(lCountEntries) + "]~Entry Format Error~" + vCatalogEntries[lCountEntries] + "\n";

									omp_set_lock(&lockList);
									ReportTimeStamp("[CatalogDatFile]", "ERROR:  _na.dat Entry [" + vCatalogEntries[lCountEntries] + "] is Not Properly Formatted");
									omp_unset_lock(&lockList);
								}
							}
						}

						//  Destroy time stamp lock
						omp_destroy_lock(&lockList);

						//  Write error file
						if (!strErrorFilePathName.empty())
						{
							//  Iterate list file vector and concatenate file
							for (long lCount = 0; lCount < vErrorEntries.size(); lCount++)
							{
								if (!vErrorEntries[lCount].empty())
									strErrorFileText += vErrorEntries[lCount];
							}

							//  Write the file
							if (!WriteFileText(strErrorFilePathName, strErrorFileText))
								ReportTimeStamp("[CatalogDatFile]", "ERROR:  Catalog Error File [" + strErrorFilePathName + "] Write Failed");
						}

						return true;
					}
					else
					{
						ReportTimeStamp("[CatalogDatFile]", "ERROR:  .dat Entries Container is Not Set");
					}
				}
				else
				{
					ReportTimeStamp("[CatalogDatFile]", "ERROR:  .dat CDS Container is Not Set");
				}
			}
			else
			{
				ReportTimeStamp("[CatalogDatFile]", "ERROR:  _na.dat Catalog File Text is Empty");
			}
		}
		else
		{
			ReportTimeStamp("[CatalogDatFile]", "ERROR:  .dat Catalog File Text is Empty");
		}
	}
	catch (exception ex)
	{
		cout << "ERROR [CatalogDatFile] Exception Code:  " << ex.what() << "\n";
	}

	return false;
}

////////////////////////////////////////////////////////////////////////////////
//
//  Catalogs a GBK file format
//
////////////////////////////////////////////////////////////////////////////////
//
//  [string&] strGBKFileText:  file text to catalog
//  [long] lIndex           :  index of the entry to set
//                         :  returns true, if successful; else, false
//         
////////////////////////////////////////////////////////////////////////////////

bool CBase123_Catalog::CatalogGBKFile(string& strGBKFileText, long lIndex)
{
	// File text lines vector
	vector<string> vLines;
	//  Entry to catalog
	CBase123_Catalog_Entry eAdd;
	//  Name/id
	string strNameID = "";
	//  The lower-case line
	string strLCLine = "";
	//  CDSs
	vector<structCDS> vCDSs;
	//  CDS name/id
	string strCDSNameID = "";
	//  In CDS string, if true
	bool bIsCDS = false;
	//  Virus, if true
	bool bIsVirus = false;

	try
	{
		//  If the gbk file text is not empty
		if (!strGBKFileText.empty())
		{
			//  Split file text
			SplitString(strGBKFileText, '\n', vLines);

			//  If lines exist
			if (vLines.size() > 0)
			{
				//  Iterate through lines, add entry to catalog
				for (int nCount = 0; nCount < vLines.size(); nCount++)
				{
					//  One-time line conversion to lower-case
					strLCLine = ConvertStringToLowerCase(vLines[nCount]);

					//  LOCUS       NC_010314               1090 bp ss-DNA     circular VRL 20-OCT-2015
					if (strLCLine.substr(0, m_strGBKFeature_Locus.length()) == m_strGBKFeature_Locus)
					{
						//  Locus fields
						vector<string> vFields;
						//  Sequence length
						long lLength = 0;

						//  Split the locus; all double spaces create empty entries, which are not pushed onto the vector
						SplitString(vLines[nCount], ' ', vFields);

						//  If the fields vector is not empty
						if (vFields.size() >= 3)
						{
							//  Accession
							eAdd.SetAccession(vFields[1]);

							//  Length
							stringstream(vFields[2]) >> lLength;
							eAdd.SetLength(lLength);
						}
						
						//  If the fields vector is not empty to 8 places, add date
						if (vFields.size() >= 8)
							eAdd.SetDate(vFields[7]);

						if (strLCLine.find(" dna ") != string::npos)
						{
							eAdd.SetMoleculeType("dna");
							eAdd.SetStrandedness("ds");
							eAdd.SetStrandednessDirection("+/-");

						}
						else if (strLCLine.find(" rna ") != string::npos)
						{
							eAdd.SetMoleculeType("rna");
							eAdd.SetStrandedness("ds");
							eAdd.SetStrandednessDirection("+/-");
						}
						else if (strLCLine.find("dna") != string::npos)
							eAdd.SetMoleculeType("dna");
						else if (strLCLine.find("rna") != string::npos)
							eAdd.SetMoleculeType("rna");

						//  Set remaining available features
						SetGBKEntryFeatures(eAdd, strLCLine);
					}

					//  Process features, definition and organism, and set name/id
					if ((strLCLine.substr(0, m_strGBKFeature_Definition.length()) == m_strGBKFeature_Definition) || (strLCLine.substr(0, m_strGBKFeature_Organism.length()) == m_strGBKFeature_Organism))
					{
						//  Next line
						string strNextLine = "";
						//  First name/id line, if true
						bool bFirstNameIDLine = true;
						//  Organism line, if true
						bool bIsOrganismLine = false;

						if (strLCLine.substr(0, m_strGBKFeature_Organism.length()) == m_strGBKFeature_Organism)
							bIsOrganismLine = true;

						//  Set initial name (all entry features have same length)
						strNameID = vLines[nCount].substr(m_strGBKFeature_Definition.length(), vLines[nCount].length() - (m_strGBKFeature_Definition.length()));

						//  Concatenate organism line(s) (definition and organism lines are the same length)
						strNextLine = vLines[nCount + 1];
						while (strNextLine.substr(0, 1) == " ")
						{
							//  If first  name/id line
							//  Increment count due to initial + 1
							if (bFirstNameIDLine)
								nCount++;

							//  Organism lines, only
							if (bIsOrganismLine)
							{
								//  If first  name/id line
								if (bFirstNameIDLine)
								{
									//  Prepend leader
									strNameID += " [";
								}
								//  Prepend leader
								else
									strNameID += " ";
							}
							//  Prepend leader
							else
								strNameID += " ";

							//  Concatenate name/id
							strNameID += strNextLine.substr(12, strNextLine.length() - 12);

							//  Increment count
							nCount++;
							//  Get next line
							strNextLine = vLines[nCount];

							//  Reset first flag
							bFirstNameIDLine = false;
						}

						//  Append closing square bracket
						if ((bIsOrganismLine) && (!bFirstNameIDLine))
							strNameID += "]";

						//  Set name/id
						if (eAdd.GetNameID().empty())
							eAdd.SetNameID(strNameID);
						else
							eAdd.SetNameID(eAdd.GetNameID() + " {" + strNameID + "}");

						//  Set remaining available features
						strNameID = ConvertStringToLowerCase(strNameID);
						SetGBKEntryFeatures(eAdd, strNameID);

						//  Set virus flag
						if (strNameID.find("viruses") != string::npos)
							bIsVirus = true;

						//  Set transcript flag
						if ((!bIsVirus) && ((strNameID.find("transcript") != string::npos) || (strNameID.find("transcribe") != string::npos) || (strNameID.find("ncrna") != string::npos) || 
							(strNameID.find("microrna") != string::npos) || (strNameID.find("non-coding rna") != string::npos) || (strNameID.find("misc_rna") != string::npos) || 
							(strNameID.find("small nucleolar rna") != string::npos) || (strNameID.find("small nuclear rna") != string::npos) || (strNameID.find("antisense rna") != string::npos) || 
							(strNameID.find("guide rna") != string::npos) || (strNameID.find("telomerase rna") != string::npos) || (strNameID.find("ribosomal rna") != string::npos) || 
							(strNameID.find("vault rna") != string::npos) || (strNameID.find(" y rna") != string::npos) || (strNameID.find("mitochondrial rna") != string::npos) || 
							(strNameID.find(" p rna") != string::npos) || (strNameID.find(" srp rna") != string::npos)))
							eAdd.SetTranscriptStatus("y");
					}
					//  Process features, version, and set versioned accession
					else if (strLCLine.substr(0, m_strGBKFeature_Version.length()) == m_strGBKFeature_Version)
					{
						//  Versioned accession
						string strVersionLine = "";
						vector<string> vVersion;

						//  Set feature, mask feature and cut final double-quote [...length() + 1]
						strVersionLine = vLines[nCount].substr(m_strGBKFeature_Version.length(), vLines[nCount].length() - (m_strGBKFeature_Version.length()));

						//  In case other features are set with this line (e.g., GI#), split and set only first entry
						SplitString(strVersionLine, ' ', vVersion);

						//  If vector contains at least one entry
						if (vVersion.size() >= 1)
							eAdd.SetAccession(ReplaceInString(vVersion[0], ".", "_", false));
					}
					//  Process features, host
					else if (strLCLine.substr(0, m_strGBKFeature_Host.length()) == m_strGBKFeature_Host)
					{
						//  Set feature, mask feature and cut final double-quote [...length() + 1]
						eAdd.SetHost(strLCLine.substr(m_strGBKFeature_Host.length(), strLCLine.length() - (m_strGBKFeature_Host.length() + 1)));
					}
					//  Process features, collection date (superior to locus date)
					else if (strLCLine.substr(0, m_strGBKFeature_CollectionDate.length()) == m_strGBKFeature_CollectionDate)
					{
						//  Set feature, mask feature and cut final double-quote [...length() + 1]
						eAdd.SetDate(strLCLine.substr(m_strGBKFeature_CollectionDate.length(), strLCLine.length() - (m_strGBKFeature_CollectionDate.length() + 1)));
					}
					//  Process features, chromosome/segment
					else if (strLCLine.substr(0, m_strGBKFeature_Chromosome.length()) == m_strGBKFeature_Chromosome)
					{
						//  Set feature, mask feature and cut final double-quote [...length() + 1] (both are same length)
						eAdd.SetChromosomeSegment(strLCLine.substr(m_strGBKFeature_Chromosome.length(), strLCLine.length() - (m_strGBKFeature_Chromosome.length() + 1)));
					}
					//  Process features, chromosome/segment
					else if (strLCLine.substr(0, m_strGBKFeature_Segment.length()) == m_strGBKFeature_Segment)
					{
						//  Set feature, mask feature and cut final double-quote [...length() + 1] (both are same length)
						eAdd.SetChromosomeSegment(strLCLine.substr(m_strGBKFeature_Segment.length(), strLCLine.length() - (m_strGBKFeature_Segment.length() + 1)));
					}
					//  Process features, country (superior to locus locale)
					else if (strLCLine.substr(0, m_strGBKFeature_Country.length()) == m_strGBKFeature_Country)
					{
						//  Set feature, mask feature and cut final double-quote [...length() + 1]
						eAdd.SetLocale(strLCLine.substr(m_strGBKFeature_Country.length(), strLCLine.length() - (m_strGBKFeature_Country.length() + 1)));
					}
					//  Process CDS features
					else if (strLCLine.substr(0, m_strGBKFeature_CDS.length()) == m_strGBKFeature_CDS)
						bIsCDS = true;

					//  CDS entry
					if (bIsCDS)
					{
						//  Process features, stop and add previous CDS, start new CDS
						if (strLCLine.substr(0, m_strGBKFeature_CDS.length()) == m_strGBKFeature_CDS)
						{
							//  Completeness
							string strCompleteness = "";
							//  Complement flag
							string strIsComplement = "n";
							//  In case of multiple starts, stops (splices/joins)
							vector<string> vCDSJoins;

							//  Set completeness
							strCompleteness = strLCLine.substr(m_strGBKFeature_CDS.length(), strLCLine.length() - m_strGBKFeature_CDS.length());

							//  Clear join characters
							if (strCompleteness.find("join") != string::npos)
							{
								strCompleteness = ReplaceInString(strCompleteness, "join(", "", true);
								strCompleteness = ReplaceInString(strCompleteness, ")", "", true);
							}

							//  Set complement flag
							if (strCompleteness.find("complement") != string::npos)
							{
								//  Set flag
								strIsComplement = "y";

								//  Clear join characters
								strCompleteness = ReplaceInString(strCompleteness, "complement(", "", true);
								strCompleteness = ReplaceInString(strCompleteness, ")", "", true);
							}

							//  Set the next CDS start/stop
							SplitString(strCompleteness, ',', vCDSJoins);

							//  Set completeness
							strCompleteness = vLines[nCount].substr(m_strGBKFeature_CDS.length(), vLines[nCount].length() - m_strGBKFeature_CDS.length());

							//  Iterate the splices/joins, if any
							for (int nCountJoins = 0; nCountJoins < vCDSJoins.size(); nCountJoins++)
							{
								//  CDS start/stop
								vector<string> vCDSStartStop;

								//  Set the next CDS start/stop
								SplitString(vCDSJoins[nCountJoins], '.', vCDSStartStop);

								//  If CDS is properly set and formatted
								if (vCDSStartStop.size() == 2)
								{
									//  CDS
									structCDS stAdd;

									stAdd.strIsComplement = strIsComplement;
									stAdd.strNameID = "";
									stringstream(vCDSStartStop[0]) >> stAdd.lStart;
									stringstream(vCDSStartStop[1]) >> stAdd.lStop;
									stAdd.strCompleteness = strCompleteness;

									//  Add CDS
									vCDSs.push_back(stAdd);
								}
								//  Clear start/stop
								vCDSStartStop.clear();
							}

							//  Clear joins
							vCDSJoins.clear();
						}
						//  Process features, CDS gene (name/id)
						else if (strLCLine.substr(0, m_strGBKFeature_Gene.length()) == m_strGBKFeature_Gene)
						{
							if (strCDSNameID.empty())
								strCDSNameID = "Gene:  " + strLCLine.substr(m_strGBKFeature_Gene.length(), strLCLine.length() - (m_strGBKFeature_Gene.length() + 1));
							else
								strCDSNameID += " [Gene:  " + strLCLine.substr(m_strGBKFeature_Gene.length(), strLCLine.length() - (m_strGBKFeature_Gene.length() + 1)) + "]";
						}
						//  Process features, CDS product (name/id)
						else if (strLCLine.substr(0, m_strGBKFeature_Product.length()) == m_strGBKFeature_Product)
						{
							//  Next line
							string strNextLine = "";

							strNextLine = ConvertStringToLowerCase(vLines[nCount + 1]);

							if (strNextLine.substr(0, m_strGBKFeature_NextFeature.length()) == m_strGBKFeature_NextFeature)
								strLCLine = strLCLine.substr(m_strGBKFeature_Product.length(), strLCLine.length() - (m_strGBKFeature_Product.length() + 1));
							else
							{
								nCount++;
								strLCLine = strLCLine.substr(m_strGBKFeature_Product.length(), strLCLine.length() - (m_strGBKFeature_Product.length()));
							}

							while (strNextLine.substr(0, m_strGBKFeature_NextFeature.length()) != m_strGBKFeature_NextFeature)
							{
								string strHoldLine = strNextLine;

								nCount++;
								strNextLine = ConvertStringToLowerCase(vLines[nCount]);

								if (strNextLine.substr(0, m_strGBKFeature_NextFeature.length()) == m_strGBKFeature_NextFeature)
									strLCLine += " " + strHoldLine.substr(21, strHoldLine.length() - 22);
								else
									strLCLine += " " + strHoldLine.substr(21, strHoldLine.length() - 21);
							}

							if (strCDSNameID.empty())
								strCDSNameID = "Product:  " + strLCLine;
							else
								strCDSNameID += " [Product:  " + strLCLine + "]";
						}
						//  Process features, CDS gene id (name/id)
						else if (strLCLine.substr(0, m_strGBKFeature_GeneID.length()) == m_strGBKFeature_GeneID)
						{
							if (strCDSNameID.empty())
								strCDSNameID = "GeneID:  " + strLCLine.substr(m_strGBKFeature_GeneID.length(), strLCLine.length() - (m_strGBKFeature_GeneID.length() + 1));
							else
								strCDSNameID += " [GeneID:  " + strLCLine.substr(m_strGBKFeature_GeneID.length(), strLCLine.length() - (m_strGBKFeature_GeneID.length() + 1)) + "]";
						}
						//  Process features, translation, add the current CDS
						else if ((strLCLine.substr(0, m_strGBKFeature_Translation.length()) == m_strGBKFeature_Translation) || (strLCLine.substr(0, m_strGBKFeature_Origin.length()) == m_strGBKFeature_Origin))
						{
							//  If the previous CDS is Not null, add it
							if (!strCDSNameID.empty())
							{
								for (int nCountCDS = 0; nCountCDS < vCDSs.size(); nCountCDS++)
								{
									if ((vCDSs[nCountCDS].lStart > 0) && (vCDSs[nCountCDS].lStop > 0) && (!vCDSs[nCountCDS].strCompleteness.empty()))
									{
										if(vCDSs.size() > 1)
											vCDSs[nCountCDS].strNameID = strCDSNameID + "{Join " + ConvertIntToString(nCountCDS + 1) + "}";
										else
											vCDSs[nCountCDS].strNameID = strCDSNameID;

										eAdd.CDSAdd(vCDSs[nCountCDS]);
									}
								}
							}

							//  Reset CDS features
							vCDSs.clear();
							strCDSNameID = "";

							//  Reset CDS flag
							bIsCDS = false;
						}
					}

					//  End catalog
					if (strLCLine.substr(0, m_strGBKFeature_Origin.length()) == m_strGBKFeature_Origin)
						break;

					strLCLine = "";
				}

				//  If virus, attempt to set viral group
				if (bIsVirus)
				{
					// i   = dsDNA (no RT)
					// ii  = ssDNA
					// iii = dsRNA
					// iv  = ssRNA+ (no RT)
					// v   = ssRNA- (no RT)
					// vi  = dsDNA (+RT)
					// vii = ssRNA (+RT)
					if ((eAdd.GetMoleculeType() == "rna") && (eAdd.GetRetroTranscriptase() == "y"))
					{
						eAdd.SetStrandedness("ss");
						eAdd.SetStrandednessDirection("+");
						eAdd.SetDNAIntermediate("y");
						eAdd.SetViralGroup("vi");
					}
					else if ((eAdd.GetMoleculeType() == "dna") && (eAdd.GetRetroTranscriptase() == "y"))
					{
						eAdd.SetStrandedness("ds");
						eAdd.SetStrandednessDirection("+/-");
						eAdd.SetRNAIntermediate("y");
						eAdd.SetViralGroup("vii");
					}
					else if ((eAdd.GetMoleculeType() == "dna") && (eAdd.GetStrandedness() == "ds"))
						eAdd.SetViralGroup("i");
					else if ((eAdd.GetMoleculeType() == "dna") && (eAdd.GetStrandedness() == "ss"))
					{
						eAdd.SetViralGroup("ii");
						eAdd.SetStrandednessDirection("+");
					}
					else if ((eAdd.GetMoleculeType() == "rna") && (eAdd.GetStrandedness() == "ds"))
						eAdd.SetViralGroup("iii");
					else if ((eAdd.GetMoleculeType() == "rna") && (eAdd.GetStrandedness() == "ss") && (eAdd.GetStrandednessDirection() == "+"))
						eAdd.SetViralGroup("iv");
					else if ((eAdd.GetMoleculeType() == "rna") && (eAdd.GetStrandedness() == "ss") && (eAdd.GetStrandednessDirection() == "-"))
						eAdd.SetViralGroup("v");

					//  In case transcript flag got set when it shouldn't (definition and organism conflict)
					eAdd.SetTranscriptStatus("");
				}

				//  Transcript
				if (eAdd.GetTranscriptStatus() == "y")
				{
					eAdd.SetMoleculeType("rna");
					eAdd.SetStrandedness("ss");
					eAdd.SetStrandednessType("l");
					eAdd.SetStrandednessDirection("-");
				}

				//  If entry has no CDS defined, the entire length is the CDS
				if (eAdd.CDSGetCount() <= 0)
					eAdd.CDSAdd("Molecular Polymer", 1, eAdd.GetLength(), "1.." + ConvertLongToString(eAdd.GetLength()));

				//  Add the entry
				SetEntryAtIndex(eAdd, lIndex);

				return true;
			}
			else
			{
				ReportTimeStamp("[CatalogGBKFile]", "ERROR:  .gbk Catalog File Container is Not Set");
			}
		}
		else
		{
			ReportTimeStamp("[CatalogGBKFile]", "ERROR:  .gbk Catalog File Text is Empty");
		}
	}
	catch (exception ex)
	{
		cout << "ERROR [CatalogGBKFile] Exception Code:  " << ex.what() << "\n";
	}

	return false;
}

////////////////////////////////////////////////////////////////////////////////
//
//  Parses a .dat CDS starts/stops string into numerics
//
////////////////////////////////////////////////////////////////////////////////
//
//  [string&] strParts          :  starts/stops string to parse
//  [vector<long>&] vStartsStops:  starts/stops vector<long> to return
//                             :  returns true, if successful; else, false
//         
////////////////////////////////////////////////////////////////////////////////

bool CBase123_Catalog::ParseStartsStops(string& strParts, vector<long>& vStartsStops)
{
	//  Numeric, string representation
	string strNumeric = "";
	//  Numeric, long representation
	long lNumeric = 0;
	//  On numeric, if true; else, on string
	bool bOnNumeric = false;

	try
	{
		//  Clear vector
		vStartsStops.clear();

		//  If parts is not empty
		if (!strParts.empty())
		{
			//  Iterate through parts string, parse out numbers by pairs
			for (int nCount = 0; nCount <= strParts.length(); nCount++)
			{
				//  On numeric
				if ((nCount < strParts.length()) && 
					((strParts[nCount] == '0') || (strParts[nCount] == '1') || (strParts[nCount] == '2') || (strParts[nCount] == '3') || (strParts[nCount] == '4') ||
					(strParts[nCount] == '5') || (strParts[nCount] == '6') || (strParts[nCount] == '7') || (strParts[nCount] == '8') || (strParts[nCount] == '9')))
				{
					bOnNumeric = true;
				}
				//  Off numeric
				else
				{
					//  If the numeric string is not empty
					if (!strNumeric.empty())
					{
						//  Convert numeric string to long
						stringstream(strNumeric) >> lNumeric;

						//  Add the numeric to the vector<long>
						vStartsStops.push_back(lNumeric);
					}

					//  Clear numeric
					strNumeric = "";
					lNumeric = 0;

					//  Off numeric
					bOnNumeric = false;
				}

				//  Concatenate numeric
				if (bOnNumeric)
					strNumeric += strParts[nCount];
			}

			return true;
		}
		else
		{
			ReportTimeStamp("[ParseStartsStops]", "ERROR:  CDS Starts/Stops is Empty");
		}
	}
	catch (exception ex)
	{
		cout << "ERROR [ParseStartsStops] Exception Code:  " << ex.what() << "\n";
	}

	return false;
}

////////////////////////////////////////////////////////////////////////////////
//
//  Gets the entries header header
//
////////////////////////////////////////////////////////////////////////////////
//
//  :  returns the entries header string
//         
////////////////////////////////////////////////////////////////////////////////

string CBase123_Catalog::GetEntryHeader()
{
	//  The header string to return
	string strHeader = "";

	try
	{
		strHeader = "Accession~";
		strHeader += "Name/ID~";
		strHeader += "Length~";
		strHeader += "Molecule Type [dna/rna]~";
		strHeader += "Strandedness [ss/ds]~";
		strHeader += "Strandedness Type [l/c]~";
		strHeader += "Strandedness Direction [+/-]~";
		strHeader += "DNA Intermediate [y/n]~";
		strHeader += "RNA Intermediate [y/n]~";
		strHeader += "Retro-Transcriptase [y/n]~";
		strHeader += "Completeness [c/nc/p]~";
		strHeader += "Date~";
		strHeader += "Host~";
		strHeader += "Host Age~";
		strHeader += "Host Gender [m/f/u]~";
		strHeader += "Chromosome/Segment~";
		strHeader += "Sero-Type~";
		strHeader += "Locale~";
		strHeader += "Satellite Status [y/n]~";
		strHeader += "Transcript Status [y/n]~";
		strHeader += "Viral Group [i-vii]";
	}
	catch (exception ex)
	{
		cout << "ERROR [GetEntryHeader] Exception Code:  " << ex.what() << "\n";
	}

	return strHeader;
}

////////////////////////////////////////////////////////////////////////////////
//
//  Gets the entry CDS header
//
////////////////////////////////////////////////////////////////////////////////
//
//  :  returns the entry CDS header
//         
////////////////////////////////////////////////////////////////////////////////

string CBase123_Catalog::GetCDSHeader()
{
	//  The entry CDS header
	string strHeader = "";
	//  The pad string
	string strPad = "";

	try
	{
		for (int nCount = 1; nCount <= 100; nCount++)
		{
			if (nCount < 10)
				strPad = "00";
			else if (nCount < 100)
				strPad = "0";

			if (nCount > 1)
				strHeader += "~";

			strHeader += "Is Complement (" + strPad + ConvertIntToString(nCount) + ")^";
			strHeader += "Name/ID (" + strPad + ConvertIntToString(nCount) + ")^";
			strHeader += "Start (" + strPad + ConvertIntToString(nCount) + ")^";
			strHeader += "Stop (" + strPad + ConvertIntToString(nCount) + ")^";
			strHeader += "Completeness (" + strPad + ConvertIntToString(nCount) + ")";
		}
	}
	catch (exception ex)
	{
		cout << "ERROR [GetCDSHeader] Exception Code:  " << ex.what() << "\n";
	}

	return strHeader;
}

////////////////////////////////////////////////////////////////////////////////
//
//  Clears catalog entries
//
////////////////////////////////////////////////////////////////////////////////
//
//  :  returns true, if successful; else, false
//         
////////////////////////////////////////////////////////////////////////////////

bool CBase123_Catalog::ClearEntries()
{
	try
	{
		m_vEntries.clear();

		return true;
	}
	catch (exception ex)
	{
		cout << "ERROR [ClearEntries] Exception Code:  " << ex.what() << "\n";
	}

	return false;
}
// Base123_Catalog_Entry.cpp

////////////////////////////////////////////////////////////////////////////////
//
//  Base123_Catalog_Entry class encapsulates all the functionality of the Base123 accession catalog
//      entry as drawn from the NCBI .gbk or .dat (influenza, only) file format(s);
//
//  Developed by Stephen Donald Huff, PhD (Stephen.Huff.3@us.af.mil)
//  Biological Informatics Group, RHDJ, 711HPW, United States Air Force Research Laboratory
//  9 July 2016
//  (All Rights Reserved)
//
////////////////////////////////////////////////////////////////////////////////

#include "F_Dist_R.h"
#include "Base123_Catalog_Entry.h"
#include "Base123_Utilities.h"

#include <sstream>

//  Initialization

////////////////////////////////////////////////////////////////////////////////
//
//  Constructs the CBase123_Catalog_Entry class object
//
////////////////////////////////////////////////////////////////////////////////
//
//
//         
////////////////////////////////////////////////////////////////////////////////

CBase123_Catalog_Entry::CBase123_Catalog_Entry()
{
	try
	{
		m_strAccession = "";
		m_strNameID = "";
		m_lLength = 0;
		m_strStrandedness = "";
		m_strStrandednessType = "";
		m_strMoleculeType = "";
		m_strStrandednessDirection = "";
		m_strDNAIntermediate = "";
		m_strRNAIntermediate = "";
		m_strRetroTranscriptase = "";
		m_strCompleteness = "";
		m_strDate = "";
		m_strHost = "";
		m_strHostAge = "";
		m_strHostGender = "";
		m_strChromosomeSegment = "";
		m_strSeroType = "";
		m_strLocale = "";
		m_strSatelliteStatus = "";
		m_strViralGroup = "";

		m_bIsSet = false;

		m_vCDSs.clear();
		m_vCDSs.reserve(100);
	}
	catch (exception ex)
	{
		cout << "ERROR [CBase123_Catalog_Entry] Exception Code:  " << ex.what() << "\n";
	}
}

////////////////////////////////////////////////////////////////////////////////
//
//  Destructs the CBase123_Catalog_Entry class object
//
////////////////////////////////////////////////////////////////////////////////
//
//
//         
////////////////////////////////////////////////////////////////////////////////

CBase123_Catalog_Entry::~CBase123_Catalog_Entry()
{
	try
	{
	}
	catch (exception ex)
	{
		cout << "ERROR [~CBase123_Catalog_Entry] Exception Code:  " << ex.what() << "\n";
	}
}

//  Interface (public)

////////////////////////////////////////////////////////////////////////////////
//
//  Gets a pointer to the CDS colletion vector<structCDS>
//
////////////////////////////////////////////////////////////////////////////////
//
//  :Returns the a pointer to the CDS collection vector<structCDS>, if successful; else, false
//         
////////////////////////////////////////////////////////////////////////////////

vector<structCDS>* CBase123_Catalog_Entry::GetCDSCollectionPointer()
{
	try
	{
		return &m_vCDSs;
	}
	catch (exception ex)
	{
		cout << "ERROR [GetCDSCollectionPointer] Exception Code:  " << ex.what() << "\n";
	}

	return NULL;
}

////////////////////////////////////////////////////////////////////////////////
//
//  Gets the CDS colletion vector<structCDS>
//
////////////////////////////////////////////////////////////////////////////////
//
//  :Returns the a pointer to the CDS collection vector<structCDS>, if successful; else, false
//         
////////////////////////////////////////////////////////////////////////////////

bool CBase123_Catalog_Entry::GetCDSCollection(vector<structCDS>& vCDSs)
{
	try
	{
		vCDSs = m_vCDSs;

		return true;
	}
	catch (exception ex)
	{
		cout << "ERROR [GetCDSCollection] Exception Code:  " << ex.what() << "\n";
	}

	return NULL;
}

////////////////////////////////////////////////////////////////////////////////
//
//  Gets the set flag
//
////////////////////////////////////////////////////////////////////////////////
//
//  :Returns the set flag
//         
////////////////////////////////////////////////////////////////////////////////

bool CBase123_Catalog_Entry::GetIsSet()
{
	try
	{
		return m_bIsSet;
	}
	catch (exception ex)
	{
		cout << "ERROR [GetIsSet] Exception Code:  " << ex.what() << "\n";
	}

	return "";
}

////////////////////////////////////////////////////////////////////////////////
//
//  Sets the set flag
//
////////////////////////////////////////////////////////////////////////////////
//
//  [bool] bIsSet:  set flag
//         
////////////////////////////////////////////////////////////////////////////////

void CBase123_Catalog_Entry::SetIsSet(bool bIsSet)
{
	try
	{
		m_bIsSet = bIsSet;
	}
	catch (exception ex)
	{
		cout << "ERROR [SetIsSet] Exception Code:  " << ex.what() << "\n";
	}

	return;
}

////////////////////////////////////////////////////////////////////////////////
//
//  Gets the accession
//
////////////////////////////////////////////////////////////////////////////////
//
//  :Returns the accession
//         
////////////////////////////////////////////////////////////////////////////////

string CBase123_Catalog_Entry::GetAccession()
{
	try
	{
		return m_strAccession;
	}
	catch (exception ex)
	{
		cout << "ERROR [GetAccession] Exception Code:  " << ex.what() << "\n";
	}

	return "";
}

////////////////////////////////////////////////////////////////////////////////
//
//  Sets the entry accession
//
////////////////////////////////////////////////////////////////////////////////
//
//  [string] strAccession:  entry accession
//         
////////////////////////////////////////////////////////////////////////////////

void CBase123_Catalog_Entry::SetAccession(string strAccession)
{
	try
	{
		m_strAccession = ScrubDelimitedEntryString(strAccession);
	}
	catch (exception ex)
	{
		cout << "ERROR [SetAccession] Exception Code:  " << ex.what() << "\n";
	}

	return;
}

////////////////////////////////////////////////////////////////////////////////
//
//  Gets the entry name/id
//
////////////////////////////////////////////////////////////////////////////////
//
//  :Returns the entry name/id
//         
////////////////////////////////////////////////////////////////////////////////

string CBase123_Catalog_Entry::GetNameID()
{
	try
	{
		return m_strNameID;
	}
	catch (exception ex)
	{
		cout << "ERROR [GetNameID] Exception Code:  " << ex.what() << "\n";
	}

	return "";
}

////////////////////////////////////////////////////////////////////////////////
//
//  Sets the entry name/id
//
////////////////////////////////////////////////////////////////////////////////
//
//  [string] strNameID:  entry name/id
//         
////////////////////////////////////////////////////////////////////////////////

void CBase123_Catalog_Entry::SetNameID(string strNameID)
{
	try
	{
		m_strNameID = ScrubDelimitedEntryString(strNameID);
	}
	catch (exception ex)
	{
		cout << "ERROR [SetNameID] Exception Code:  " << ex.what() << "\n";
	}

	return;
}

////////////////////////////////////////////////////////////////////////////////
//
//  Gets the entry  length
//
////////////////////////////////////////////////////////////////////////////////
//
//  :Returns the entry length
//         
////////////////////////////////////////////////////////////////////////////////

long CBase123_Catalog_Entry::GetLength()
{
	try
	{
		return m_lLength;
	}
	catch (exception ex)
	{
		cout << "ERROR [GetLength] Exception Code:  " << ex.what() << "\n";
	}

	return 0;
}

////////////////////////////////////////////////////////////////////////////////
//
//  Sets the entry length
//
////////////////////////////////////////////////////////////////////////////////
//
//  [long] lLength:  entry length
//         
////////////////////////////////////////////////////////////////////////////////

void CBase123_Catalog_Entry::SetLength(long lLength)
{
	try
	{
		m_lLength = lLength;
	}
	catch (exception ex)
	{
		cout << "ERROR [SetLength] Exception Code:  " << ex.what() << "\n";
	}

	return;
}

////////////////////////////////////////////////////////////////////////////////
//
//  Gets the entry molecule type [dna/rna]
//
////////////////////////////////////////////////////////////////////////////////
//
//  :Returns the entry molecule type [dna/rna]
//         
////////////////////////////////////////////////////////////////////////////////

string CBase123_Catalog_Entry::GetMoleculeType()
{
	try
	{
		return m_strMoleculeType;
	}
	catch (exception ex)
	{
		cout << "ERROR [GetMoleculeType] Exception Code:  " << ex.what() << "\n";
	}

	return "";
}

////////////////////////////////////////////////////////////////////////////////
//
//  Sets the entry molecule type [dna/rna[
//
////////////////////////////////////////////////////////////////////////////////
//
//  [string] strMoleculeType:  entry molecule type [dna/rna]
//         
////////////////////////////////////////////////////////////////////////////////

void CBase123_Catalog_Entry::SetMoleculeType(string strMoleculeType)
{
	try
	{
		strMoleculeType = ConvertStringToLowerCase(strMoleculeType);

		if((strMoleculeType.empty()) || (strMoleculeType == "dna") || (strMoleculeType == "rna"))
			m_strMoleculeType = strMoleculeType;
		else
			ReportTimeStamp("[SetMoleculeType]", "ERROR:  Molecule type [" + strMoleculeType + "] must be either [dna] or [rna]");
	}
	catch (exception ex)
	{
		cout << "ERROR [SetMoleculeType] Exception Code:  " << ex.what() << "\n";
	}

	return;
}

////////////////////////////////////////////////////////////////////////////////
//
//  Gets the entry strandedness [ds/ss]
//
////////////////////////////////////////////////////////////////////////////////
//
//  :Returns the entry strandedness [ds/ss]
//         
////////////////////////////////////////////////////////////////////////////////

string CBase123_Catalog_Entry::GetStrandedness()
{
	try
	{
		return m_strStrandedness;
	}
	catch (exception ex)
	{
		cout << "ERROR [GetStrandedness] Exception Code:  " << ex.what() << "\n";
	}

	return "";
}

////////////////////////////////////////////////////////////////////////////////
//
//  Sets the entry strandedness [ds/ss]
//
////////////////////////////////////////////////////////////////////////////////
//
//  [strStrandedness] strDivison:  entry strandedness [ds/ss]
//         
////////////////////////////////////////////////////////////////////////////////

void CBase123_Catalog_Entry::SetStrandedness(string strStrandedness)
{
	try
	{
		strStrandedness = ConvertStringToLowerCase(strStrandedness);

		if((strStrandedness.empty()) || (strStrandedness == "ds") || (strStrandedness == "ss"))
			m_strStrandedness = strStrandedness;
		else
			ReportTimeStamp("[SetStrandedness]", "ERROR:  Strandedness [" + strStrandedness + "] must be either double-stranded [ds] or single-stranded [ss]");
	}
	catch (exception ex)
	{
		cout << "ERROR [SetStrandedness] Exception Code:  " << ex.what() << "\n";
	}

	return;
}

////////////////////////////////////////////////////////////////////////////////
//
//  Gets the entry strandedness type [c/l]
//
////////////////////////////////////////////////////////////////////////////////
//
//  :Returns the entry strandedness type [c/l]
//         
////////////////////////////////////////////////////////////////////////////////

string CBase123_Catalog_Entry::GetStrandednessType()
{
	try
	{
		return m_strStrandednessType;
	}
	catch (exception ex)
	{
		cout << "ERROR [GetStrandednessType] Exception Code:  " << ex.what() << "\n";
	}

	return "";
}

////////////////////////////////////////////////////////////////////////////////
//
//  Sets the entry strandedness type [c/l]
//
////////////////////////////////////////////////////////////////////////////////
//
//  [strStrandednessType] strDivison:  entry strandedness type [c/l]
//         
////////////////////////////////////////////////////////////////////////////////

void CBase123_Catalog_Entry::SetStrandednessType(string strStrandednessType)
{
	try
	{
		strStrandednessType = ConvertStringToLowerCase(strStrandednessType);

		if ((strStrandednessType.empty()) || (strStrandednessType == "c") || (strStrandednessType == "l"))
			m_strStrandednessType = strStrandednessType;
		else
			ReportTimeStamp("[SetStrandednessType]", "ERROR:  Strandedness-type [" + strStrandednessType + "] must be either circular [c] or linear [l]");
	}
	catch (exception ex)
	{
		cout << "ERROR [SetStrandednessType] Exception Code:  " << ex.what() << "\n";
	}

	return;
}

////////////////////////////////////////////////////////////////////////////////
//
//  Gets the entry strandedness direction [+/-]
//
////////////////////////////////////////////////////////////////////////////////
//
//  :Returns the entry strandedness direction [+/-]
//         
////////////////////////////////////////////////////////////////////////////////

string CBase123_Catalog_Entry::GetStrandednessDirection()
{
	try
	{
		return m_strStrandednessDirection;
	}
	catch (exception ex)
	{
		cout << "ERROR [GetStrandednessDirection] Exception Code:  " << ex.what() << "\n";
	}

	return "";
}

////////////////////////////////////////////////////////////////////////////////
//
//  Sets the entry strandedness direction [+/-]
//
////////////////////////////////////////////////////////////////////////////////
//
//  [string] strStrandednessDirection:  entry strandedness direction [+/-]
//         
////////////////////////////////////////////////////////////////////////////////

void CBase123_Catalog_Entry::SetStrandednessDirection(string strStrandednessDirection)
{
	try
	{
		strStrandednessDirection = ConvertStringToLowerCase(strStrandednessDirection);

		if ((strStrandednessDirection.empty()) || (strStrandednessDirection == "+") || (strStrandednessDirection == "-") || (strStrandednessDirection == "+/-"))
			m_strStrandednessDirection = strStrandednessDirection;
		else
			ReportTimeStamp("[SetStrandednessDirection]", "ERROR:  Strandedness-direction [" + strStrandednessDirection + "] must be either positive [+] or negative [-] or both [+/-]");
	}
	catch (exception ex)
	{
		cout << "ERROR [SetStrandednessDirection] Exception Code:  " << ex.what() << "\n";
	}

	return;
}

////////////////////////////////////////////////////////////////////////////////
//
//  Gets the entry DNA intermediate (viruses) [y/n]
//
////////////////////////////////////////////////////////////////////////////////
//
//  :Returns the entry DNA intermediate (viruses) [y/n]
//         
////////////////////////////////////////////////////////////////////////////////

string CBase123_Catalog_Entry::GetDNAIntermediate()
{
	try
	{
		return m_strDNAIntermediate;
	}
	catch (exception ex)
	{
		cout << "ERROR [GetDNAIntermediate] Exception Code:  " << ex.what() << "\n";
	}

	return "";
}

////////////////////////////////////////////////////////////////////////////////
//
//  Sets the entry DNA intermediate (viruses) [y/n]
//
////////////////////////////////////////////////////////////////////////////////
//
//  [string] strDNAIntermediate:  entry DNA intermediate [y/n]
//         
////////////////////////////////////////////////////////////////////////////////

void CBase123_Catalog_Entry::SetDNAIntermediate(string strDNAIntermediate)
{
	try
	{
		strDNAIntermediate = ConvertStringToLowerCase(strDNAIntermediate);

		if ((strDNAIntermediate.empty()) || (strDNAIntermediate == "y") || (strDNAIntermediate == "n"))
			m_strDNAIntermediate = strDNAIntermediate;
		else
			ReportTimeStamp("[SetDNAIntermediate]", "ERROR:  DNA-intermediate [" + strDNAIntermediate + "] must be either yes [y] or no [n]");
	}
	catch (exception ex)
	{
		cout << "ERROR [SetDNAIntermediate] Exception Code:  " << ex.what() << "\n";
	}

	return;
}

////////////////////////////////////////////////////////////////////////////////
//
//  Gets the entry RNA intermediate (viruses) [y/n]
//
////////////////////////////////////////////////////////////////////////////////
//
//  :Returns the entry RNA intermediate (viruses) [y/n]
//         
////////////////////////////////////////////////////////////////////////////////

string CBase123_Catalog_Entry::GetRNAIntermediate()
{
	try
	{
		return m_strRNAIntermediate;
	}
	catch (exception ex)
	{
		cout << "ERROR [GetRNAIntermediate] Exception Code:  " << ex.what() << "\n";
	}

	return "";
}

////////////////////////////////////////////////////////////////////////////////
//
//  Sets the entry RNA intermediate (viruses) [y/n]
//
////////////////////////////////////////////////////////////////////////////////
//
//  [string] strRNAIntermediate:  entry RNA intermediate [y/n]
//         
////////////////////////////////////////////////////////////////////////////////

void CBase123_Catalog_Entry::SetRNAIntermediate(string strRNAIntermediate)
{
	try
	{
		strRNAIntermediate = ConvertStringToLowerCase(strRNAIntermediate);

		if ((strRNAIntermediate.empty()) || (strRNAIntermediate == "y") || (strRNAIntermediate == "n"))
			m_strRNAIntermediate = strRNAIntermediate;
		else
			ReportTimeStamp("[SetRNAIntermediate]", "ERROR:  RNA-intermediate [" + strRNAIntermediate + "] must be either yes [y] or no [n]");
	}
	catch (exception ex)
	{
		cout << "ERROR [SetRNAIntermediate] Exception Code:  " << ex.what() << "\n";
	}

	return;
}

////////////////////////////////////////////////////////////////////////////////
//
//  Gets the entry retro-transcriptase (viruses) [y/n]
//
////////////////////////////////////////////////////////////////////////////////
//
//  :Returns the entry retro-transcriptase (viruses) [y/n]
//         
////////////////////////////////////////////////////////////////////////////////

string CBase123_Catalog_Entry::GetRetroTranscriptase()
{
	try
	{
		return m_strRetroTranscriptase;
	}
	catch (exception ex)
	{
		cout << "ERROR [GetRetroTranscriptase] Exception Code:  " << ex.what() << "\n";
	}

	return "";
}

////////////////////////////////////////////////////////////////////////////////
//
//  Sets the entry retro-transcriptase (viruses) [y/n]
//
////////////////////////////////////////////////////////////////////////////////
//
//  [string] strRetroTranscriptase:  entry retro-transcriptase [y/n]
//         
////////////////////////////////////////////////////////////////////////////////

void CBase123_Catalog_Entry::SetRetroTranscriptase(string strRetroTranscriptase)
{
	try
	{
		strRetroTranscriptase = ConvertStringToLowerCase(strRetroTranscriptase);

		if ((strRetroTranscriptase.empty()) || (strRetroTranscriptase == "y") || (strRetroTranscriptase == "n"))
			m_strRetroTranscriptase = strRetroTranscriptase;
		else
			ReportTimeStamp("[SetRetroTranscriptase]", "ERROR:  Retro-transcriptase [" + strRetroTranscriptase + "] must be either yes [y] or no [n]");
	}
	catch (exception ex)
	{
		cout << "ERROR [SetRetroTranscriptase] Exception Code:  " << ex.what() << "\n";
	}

	return;
}

////////////////////////////////////////////////////////////////////////////////
//
//  Gets the entry completeness [c/p/nc]
//
////////////////////////////////////////////////////////////////////////////////
//
//  :Returns the entry completeness [c/p/nc]
//         
////////////////////////////////////////////////////////////////////////////////

string CBase123_Catalog_Entry::GetCompleteness()
{
	try
	{
		return m_strCompleteness;
	}
	catch (exception ex)
	{
		cout << "ERROR [GetCompleteness] Exception Code:  " << ex.what() << "\n";
	}

	return "";
}

////////////////////////////////////////////////////////////////////////////////
//
//  Sets the entry completeness [c/p/nc]
//
////////////////////////////////////////////////////////////////////////////////
//
//  [string] strCompleteness:  entry completeness [c/p/nc]
//         
////////////////////////////////////////////////////////////////////////////////

void CBase123_Catalog_Entry::SetCompleteness(string strCompleteness)
{
	try
	{
		strCompleteness = ConvertStringToLowerCase(strCompleteness);

		if ((strCompleteness.empty()) || (strCompleteness == "c") || (strCompleteness == "p") || (strCompleteness == "nc"))
			m_strCompleteness = strCompleteness;
		else
			ReportTimeStamp("[SetCompleteness]", "ERROR:  Completeness [" + strCompleteness + "] must be either complete [c] or partial [p] nearly-complete [nc]");
	}
	catch (exception ex)
	{
		cout << "ERROR [SetCompleteness] Exception Code:  " << ex.what() << "\n";
	}

	return;
}

////////////////////////////////////////////////////////////////////////////////
//
//  Gets the entry date
//
////////////////////////////////////////////////////////////////////////////////
//
//  :Returns the entry date
//         
////////////////////////////////////////////////////////////////////////////////

string CBase123_Catalog_Entry::GetDate()
{
	try
	{
		return m_strDate;
	}
	catch (exception ex)
	{
		cout << "ERROR [GetDate] Exception Code:  " << ex.what() << "\n";
	}

	return "";
}

////////////////////////////////////////////////////////////////////////////////
//
//  Sets the entry date
//
////////////////////////////////////////////////////////////////////////////////
//
//  [string] strDate:  entry date
//         
////////////////////////////////////////////////////////////////////////////////

void CBase123_Catalog_Entry::SetDate(string strDate)
{
	try
	{
		m_strDate = ScrubDelimitedEntryString(strDate);
	}
	catch (exception ex)
	{
		cout << "ERROR [SetDate] Exception Code:  " << ex.what() << "\n";
	}

	return;
}

////////////////////////////////////////////////////////////////////////////////
//
//  Gets the entry host
//
////////////////////////////////////////////////////////////////////////////////
//
//  :Returns the entry host
//         
////////////////////////////////////////////////////////////////////////////////

string CBase123_Catalog_Entry::GetHost()
{
	try
	{
		return m_strHost;
	}
	catch (exception ex)
	{
		cout << "ERROR [GetHost] Exception Code:  " << ex.what() << "\n";
	}

	return "";
}

////////////////////////////////////////////////////////////////////////////////
//
//  Sets the entry host
//
////////////////////////////////////////////////////////////////////////////////
//
//  [string] strHost:  entry host
//         
////////////////////////////////////////////////////////////////////////////////

void CBase123_Catalog_Entry::SetHost(string strHost)
{
	try
	{
		m_strHost = ScrubDelimitedEntryString(strHost);
	}
	catch (exception ex)
	{
		cout << "ERROR [SetHost] Exception Code:  " << ex.what() << "\n";
	}

	return;
}

////////////////////////////////////////////////////////////////////////////////
//
//  Gets the entry host age
//
////////////////////////////////////////////////////////////////////////////////
//
//  :Returns the entry host age
//         
////////////////////////////////////////////////////////////////////////////////

string CBase123_Catalog_Entry::GetHostAge()
{
	try
	{
		return m_strHostAge;
	}
	catch (exception ex)
	{
		cout << "ERROR [GetHostAge] Exception Code:  " << ex.what() << "\n";
	}

	return "";
}

////////////////////////////////////////////////////////////////////////////////
//
//  Sets the entry host age
//
////////////////////////////////////////////////////////////////////////////////
//
//  [string] strHostAge:  entry host age
//         
////////////////////////////////////////////////////////////////////////////////

void CBase123_Catalog_Entry::SetHostAge(string strHostAge)
{
	try
	{
		m_strHostAge = ScrubDelimitedEntryString(strHostAge);
	}
	catch (exception ex)
	{
		cout << "ERROR [SetHostAge] Exception Code:  " << ex.what() << "\n";
	}

	return;
}

////////////////////////////////////////////////////////////////////////////////
//
//  Gets the entry host gender [m/f/u]
//
////////////////////////////////////////////////////////////////////////////////
//
//  :Returns the entry host gender [m/f/u]
//         
////////////////////////////////////////////////////////////////////////////////

string CBase123_Catalog_Entry::GetHostGender()
{
	try
	{
		return m_strHostGender;
	}
	catch (exception ex)
	{
		cout << "ERROR [GetHostGender] Exception Code:  " << ex.what() << "\n";
	}

	return "";
}

////////////////////////////////////////////////////////////////////////////////
//
//  Sets the entry host gender [m/f/u]
//
////////////////////////////////////////////////////////////////////////////////
//
//  [string] strHostGender:  entry host gender [m/f/u]
//         
////////////////////////////////////////////////////////////////////////////////

void CBase123_Catalog_Entry::SetHostGender(string strHostGender)
{
	try
	{
		strHostGender = ConvertStringToLowerCase(strHostGender);

		if ((strHostGender.empty()) || (strHostGender == "m") || (strHostGender == "f") || (strHostGender == "u"))
			m_strHostGender = strHostGender;
		else
			ReportTimeStamp("[SetHostGender]", "ERROR:  Host gender [" + strHostGender + "] must be either male [m], female [f], or unknown [u]");
	}
	catch (exception ex)
	{
		cout << "ERROR [SetHostGender] Exception Code:  " << ex.what() << "\n";
	}

	return;
}

////////////////////////////////////////////////////////////////////////////////
//
//  Gets the entry segment
//
////////////////////////////////////////////////////////////////////////////////
//
//  :Returns the entry segment
//         
////////////////////////////////////////////////////////////////////////////////

string CBase123_Catalog_Entry::GetChromosomeSegment()
{
	try
	{
		return m_strChromosomeSegment;
	}
	catch (exception ex)
	{
		cout << "ERROR [GetChromosomeSegment] Exception Code:  " << ex.what() << "\n";
	}

	return 0;
}

////////////////////////////////////////////////////////////////////////////////
//
//  Sets the entry segment
//
////////////////////////////////////////////////////////////////////////////////
//
//  [int] nChromosomeSegment:  entry segment
//         
////////////////////////////////////////////////////////////////////////////////

void CBase123_Catalog_Entry::SetChromosomeSegment(string strChromosomeSegment)
{
	try
	{
		m_strChromosomeSegment = ScrubDelimitedEntryString(strChromosomeSegment);
	}
	catch (exception ex)
	{
		cout << "ERROR [SetChromosomeSegment] Exception Code:  " << ex.what() << "\n";
	}

	return;
}

////////////////////////////////////////////////////////////////////////////////
//
//  Gets the entry sero-type
//
////////////////////////////////////////////////////////////////////////////////
//
//  :Returns the entry sero-type
//         
////////////////////////////////////////////////////////////////////////////////

string CBase123_Catalog_Entry::GetSeroType()
{
	try
	{
		return m_strSeroType;
	}
	catch (exception ex)
	{
		cout << "ERROR [GetSeroType] Exception SeroType:  " << ex.what() << "\n";
	}

	return "";
}

////////////////////////////////////////////////////////////////////////////////
//
//  Sets the entry sero-type
//
////////////////////////////////////////////////////////////////////////////////
//
//  [string] strSeroType:  entry sero-type
//         
////////////////////////////////////////////////////////////////////////////////

void CBase123_Catalog_Entry::SetSeroType(string strSeroType)
{
	try
	{
		m_strSeroType = ScrubDelimitedEntryString(strSeroType);
	}
	catch (exception ex)
	{
		cout << "ERROR [SetSeroType] Exception Code:  " << ex.what() << "\n";
	}

	return;
}

////////////////////////////////////////////////////////////////////////////////
//
//  Gets the entry locale
//
////////////////////////////////////////////////////////////////////////////////
//
//  :Returns the entry locale
//         
////////////////////////////////////////////////////////////////////////////////

string CBase123_Catalog_Entry::GetLocale()
{
	try
	{
		return m_strLocale;
	}
	catch (exception ex)
	{
		cout << "ERROR [GetLocale] Exception Code:  " << ex.what() << "\n";
	}

	return "";
}

////////////////////////////////////////////////////////////////////////////////
//
//  Sets the entry locale
//
////////////////////////////////////////////////////////////////////////////////
//
//  [string] strLocale:  entry locale
//         
////////////////////////////////////////////////////////////////////////////////

void CBase123_Catalog_Entry::SetLocale(string strLocale)
{
	try
	{
		m_strLocale = ScrubDelimitedEntryString(strLocale);
	}
	catch (exception ex)
	{
		cout << "ERROR [SetLocale] Exception Code:  " << ex.what() << "\n";
	}

	return;
}

////////////////////////////////////////////////////////////////////////////////
//
//  Gets the entry satellite status (viruses) [y/n]
//
////////////////////////////////////////////////////////////////////////////////
//
//  :Returns the entry satellite status (viruses) [y/n]
//         
////////////////////////////////////////////////////////////////////////////////

string CBase123_Catalog_Entry::GetSatelliteStatus()
{
	try
	{
		return m_strSatelliteStatus;
	}
	catch (exception ex)
	{
		cout << "ERROR [GetSatelliteStatus] Exception Code:  " << ex.what() << "\n";
	}

	return "";
}

////////////////////////////////////////////////////////////////////////////////
//
//  Sets the entry satellite status (viruses) [y/n]
//
////////////////////////////////////////////////////////////////////////////////
//
//  [string] strSatelliteStatus:  entry satellite status [y/n]
//         
////////////////////////////////////////////////////////////////////////////////

void CBase123_Catalog_Entry::SetSatelliteStatus(string strSatelliteStatus)
{
	try
	{
		strSatelliteStatus = ConvertStringToLowerCase(strSatelliteStatus);

		if ((strSatelliteStatus.empty()) || (strSatelliteStatus == "y") || (strSatelliteStatus == "n"))
			m_strSatelliteStatus = strSatelliteStatus;
		else
			ReportTimeStamp("[SetSatelliteStatus]", "ERROR:  DNA-intermediate [" + strSatelliteStatus + "] must be either yes [y] or no [n]");
	}
	catch (exception ex)
	{
		cout << "ERROR [SetSatelliteStatus] Exception Code:  " << ex.what() << "\n";
	}

	return;
}

////////////////////////////////////////////////////////////////////////////////
//
//  Gets the entry transcript status (viruses) [y/n]
//
////////////////////////////////////////////////////////////////////////////////
//
//  :Returns the entry transcript status (viruses) [y/n]
//         
////////////////////////////////////////////////////////////////////////////////

string CBase123_Catalog_Entry::GetTranscriptStatus()
{
	try
	{
		return m_strTranscriptStatus;
	}
	catch (exception ex)
	{
		cout << "ERROR [GetTranscriptStatus] Exception Code:  " << ex.what() << "\n";
	}

	return "";
}

////////////////////////////////////////////////////////////////////////////////
//
//  Sets the entry transcript status (viruses) [y/n]
//
////////////////////////////////////////////////////////////////////////////////
//
//  [string] strTranscriptStatus:  entry transcript status [y/n]
//         
////////////////////////////////////////////////////////////////////////////////

void CBase123_Catalog_Entry::SetTranscriptStatus(string strTranscriptStatus)
{
	try
	{
		strTranscriptStatus = ConvertStringToLowerCase(strTranscriptStatus);

		if ((strTranscriptStatus.empty()) || (strTranscriptStatus == "y") || (strTranscriptStatus == "n"))
			m_strTranscriptStatus = strTranscriptStatus;
		else
			ReportTimeStamp("[SetTranscriptStatus]", "ERROR:  DNA-intermediate [" + strTranscriptStatus + "] must be either yes [y] or no [n]");
	}
	catch (exception ex)
	{
		cout << "ERROR [SetTranscriptStatus] Exception Code:  " << ex.what() << "\n";
	}

	return;
}

////////////////////////////////////////////////////////////////////////////////
//
//  Gets the viral group
//
////////////////////////////////////////////////////////////////////////////////
//
//  :Returns the viral group
//         
////////////////////////////////////////////////////////////////////////////////

string CBase123_Catalog_Entry::GetViralGroup()
{
	try
	{
		return m_strViralGroup;
	}
	catch (exception ex)
	{
		cout << "ERROR [GetViralGroup] Exception Code:  " << ex.what() << "\n";
	}

	return "";
}

////////////////////////////////////////////////////////////////////////////////
//
//  Sets the viral group
//
////////////////////////////////////////////////////////////////////////////////
//
//  [string] strViralGroup:  viral group
//         
////////////////////////////////////////////////////////////////////////////////

void CBase123_Catalog_Entry::SetViralGroup(string strViralGroup)
{
	try
	{
		strViralGroup = ConvertStringToLowerCase(strViralGroup);

		if ((strViralGroup.empty()) || (strViralGroup == "i") || (strViralGroup == "ii") || (strViralGroup == "iii") || (strViralGroup == "iv") || 
			(strViralGroup == "v") || (strViralGroup == "vi") || (strViralGroup == "vii"))
			m_strViralGroup = strViralGroup;
		else
			ReportTimeStamp("[SetHostGender]", "ERROR:  Viral group [" + strViralGroup + "] must be [i, ii, iii, iv, v, vi, vii]");
	}
	catch (exception ex)
	{
		cout << "ERROR [SetViralGroup] Exception Code:  " << ex.what() << "\n";
	}

	return;
}

////////////////////////////////////////////////////////////////////////////////
//
//  Adds a CDS to the CDS list
//
////////////////////////////////////////////////////////////////////////////////
//
//  [structCDS] structCDS:  CDS structure to add
//         
////////////////////////////////////////////////////////////////////////////////

void CBase123_Catalog_Entry::CDSAdd(structCDS structCDS)
{
	try
	{
		structCDS.strNameID = ScrubDelimitedEntryString(structCDS.strNameID);
		structCDS.strCompleteness = ScrubDelimitedEntryString(structCDS.strCompleteness);

		m_vCDSs.push_back(structCDS);
	}
	catch (exception ex)
	{
		cout << "ERROR [CDSAdd] Exception Code:  " << ex.what() << "\n";
	}

	return;
}

////////////////////////////////////////////////////////////////////////////////
//
//  Adds a CDS to the CDS list
//
////////////////////////////////////////////////////////////////////////////////
//
//  [string] strNameID      :  CDS name/id
//  [long] lStart           :  CDS start
//  [long] lStop            :  CDS stop
//  [string] strCompleteness:  CDS completeness
//         
////////////////////////////////////////////////////////////////////////////////

void CBase123_Catalog_Entry::CDSAdd(string strNameID, long lStart, long lStop, string strCompleteness)
{
	try
	{
		structCDS structCDS;

		structCDS.strNameID = ScrubDelimitedEntryString(strNameID);
		structCDS.lStart = lStart;
		structCDS.lStop = lStop;
		structCDS.strCompleteness = ScrubDelimitedEntryString(strCompleteness);

		m_vCDSs.push_back(structCDS);
	}
	catch (exception ex)
	{
		cout << "ERROR [CDSAdd] Exception Code:  " << ex.what() << "\n";
	}

	return;
}

////////////////////////////////////////////////////////////////////////////////
//
//  Removes the CDS structure at the indicated index
//
////////////////////////////////////////////////////////////////////////////////
//
//  [int] nIndex:  index of the CDS to modify
//         
////////////////////////////////////////////////////////////////////////////////

void CBase123_Catalog_Entry::CDSRemove(int nIndex)
{
	try
	{
		//  If the index is in range
		if ((nIndex >= 0) && (nIndex < m_vCDSs.size()))
		{
			m_vCDSs.erase(m_vCDSs.begin() + nIndex);
		}
		else
		{
			ReportTimeStamp("[CDSRemove]", "ERROR:  CDS Index [" + ConvertIntToString(nIndex) + "] is Out of Range [0:" + ConvertIntToString((int)m_vCDSs.size() - 1) + "]");
		}
	}
	catch (exception ex)
	{
		cout << "ERROR [CDSRemove] Exception Code:  " << ex.what() << "\n";
	}

	return;
}

////////////////////////////////////////////////////////////////////////////////
//
//  Gets the CDS structure at the indicated index
//
////////////////////////////////////////////////////////////////////////////////
//
//  [int] nIndex  :  index of the CDS to get
//  [structCDS&] stGet:  CDS structure returned
//
////////////////////////////////////////////////////////////////////////////////

 void CBase123_Catalog_Entry::CDSGet(int nIndex, structCDS& stGet)
{
	try
	{
		//  If the index is in range
		if ((nIndex >= 0) && (nIndex < m_vCDSs.size()))
		{
			stGet = m_vCDSs[nIndex];
		}
		else
		{
			ReportTimeStamp("[CDSGet]", "ERROR:  CDS Index [" + ConvertIntToString(nIndex) + "] is Out of Range [0:" + ConvertIntToString((int)m_vCDSs.size() - 1) + "]");
		}
	}
	catch (exception ex)
	{
		cout << "ERROR [CDSGet] Exception Code:  " << ex.what() << "\n";
	}

	return;
}

////////////////////////////////////////////////////////////////////////////////
//
//  Gets the CDS structure at the indicated index
//
////////////////////////////////////////////////////////////////////////////////
//
//  [int] nIndex:  index of the CDS to get
//             :returns the start/stop structure at the indicated index
//         
////////////////////////////////////////////////////////////////////////////////

 structCDS* CBase123_Catalog_Entry::CDSGet(int nIndex)
 {
	 try
	 {
		 //  If the index is in range
		 if ((nIndex >= 0) && (nIndex < m_vCDSs.size()))
			 return &m_vCDSs[nIndex];
		 else
		 {
			 ReportTimeStamp("[CDSGet]", "ERROR:  CDS Index [" + ConvertIntToString(nIndex) + "] is Out of Range [0:" + ConvertIntToString((int)m_vCDSs.size() - 1) + "]");
		 }
	 }
	 catch (exception ex)
	 {
		 cout << "ERROR [CDSGet] Exception Code:  " << ex.what() << "\n";
	 }

	 return NULL;
 }

////////////////////////////////////////////////////////////////////////////////
//
//  Gets the CDS collection count within the entry
//
////////////////////////////////////////////////////////////////////////////////
//
//  :returns the CDS collection count
//         
////////////////////////////////////////////////////////////////////////////////

 int CBase123_Catalog_Entry::CDSGetCount()
 {
	 try
	 {
		 return (int) m_vCDSs.size();
	 }
	 catch (exception ex)
	 {
		 cout << "ERROR [CDSGetCount] Exception Code:  " << ex.what() << "\n";
	 }

	 return 0;
 }

////////////////////////////////////////////////////////////////////////////////
//
//  Gets the complement flag of the CDS at the indicated index
//
////////////////////////////////////////////////////////////////////////////////
//
//  [int] nIndex:  index of the CDS to get
//             :  returns the complement flag at the indicated index
//         
////////////////////////////////////////////////////////////////////////////////

 string CBase123_Catalog_Entry::CDSGetIsComplement(int nIndex)
 {
	 try
	 {
		 //  If the index is in range
		 if ((nIndex >= 0) && (nIndex < m_vCDSs.size()))
		 {
			 return m_vCDSs[nIndex].strIsComplement;
		 }
		 else
		 {
			 ReportTimeStamp("[CDSGetIsComplement]", "ERROR:  CDS Index [" + ConvertIntToString(nIndex) + "] is Out of Range [0:" + ConvertIntToString((int)m_vCDSs.size() - 1) + "]");
		 }
	 }
	 catch (exception ex)
	 {
		 cout << "ERROR [CDSGetIsComplement] Exception Code:  " << ex.what() << "\n";
	 }

	 return "";
 }

////////////////////////////////////////////////////////////////////////////////
//
//  Sets the complement flag of the CDS at the indicated index
//
////////////////////////////////////////////////////////////////////////////////
//
//  [int] nIndex      :  index of the CDS to modify
//  [string] strIsComplement:  CDS complement flag to set
//         
////////////////////////////////////////////////////////////////////////////////

 void CBase123_Catalog_Entry::CDSSetIsComplement(int nIndex, string strIsComplement)
 {
	 try
	 {
		 //  If the index is in range
		 if ((nIndex >= 0) && (nIndex < m_vCDSs.size()))
			 m_vCDSs[nIndex].strIsComplement = strIsComplement;
		 else
		 {
			 ReportTimeStamp("[CDSSetIsComplement]", "ERROR:  CDS Index [" + ConvertIntToString(nIndex) + "] is Out of Range [0:" + ConvertIntToString((int)m_vCDSs.size() - 1) + "]");
		 }
	 }
	 catch (exception ex)
	 {
		 cout << "ERROR [CDSSetIsComplement] Exception Code:  " << ex.what() << "\n";
	 }
 }

////////////////////////////////////////////////////////////////////////////////
//
//  Gets the name/id of the CDS at the indicated index
//
////////////////////////////////////////////////////////////////////////////////
//
//  [int] nIndex:  index of the CDS to get
//             :  returns the name/id at the indicated index
//         
////////////////////////////////////////////////////////////////////////////////

string CBase123_Catalog_Entry::CDSGetNameID(int nIndex)
{
	try
	{
		//  If the index is in range
		if ((nIndex >= 0) && (nIndex < m_vCDSs.size()))
			return m_vCDSs[nIndex].strNameID;
		else
		{
			ReportTimeStamp("[CDSGetNameID]", "ERROR:  CDS Index [" + ConvertIntToString(nIndex) + "] is Out of Range [0:" + ConvertIntToString((int)m_vCDSs.size() - 1) + "]");
		}
	}
	catch (exception ex)
	{
		cout << "ERROR [CDSGetNameID] Exception Code:  " << ex.what() << "\n";
	}

	return "";
}

////////////////////////////////////////////////////////////////////////////////
//
//  Sets the name/id of the CDS at the indicated index
//
////////////////////////////////////////////////////////////////////////////////
//
//  [int] nIndex      :  index of the CDS to modify
//  [string] strNameID:  CDS name/id to set
//         
////////////////////////////////////////////////////////////////////////////////

void CBase123_Catalog_Entry::CDSSetNameID(int nIndex, string strNameID)
{
	try
	{
		//  If the index is in range
		if ((nIndex >= 0) && (nIndex < m_vCDSs.size()))
			m_vCDSs[nIndex].strNameID = strNameID;
		else
		{
			ReportTimeStamp("[CDSSetNameID]", "ERROR:  CDS Index [" + ConvertIntToString(nIndex) + "] is Out of Range [0:" + ConvertIntToString((int)m_vCDSs.size() - 1) + "]");
		}
	}
	catch (exception ex)
	{
		cout << "ERROR [CDSSetNameID] Exception Code:  " << ex.what() << "\n";
	}
}

////////////////////////////////////////////////////////////////////////////////
//
//  Gets the start of the CDS at the indicated index
//
////////////////////////////////////////////////////////////////////////////////
//
//  [int] nIndex:  index of the CDS to get
//             :  returns the saart at the indicated index
//         
////////////////////////////////////////////////////////////////////////////////

long CBase123_Catalog_Entry::CDSGetStart(int nIndex)
{
	try
	{
		//  If the index is in range
		if ((nIndex >= 0) && (nIndex < m_vCDSs.size()))
			return m_vCDSs[nIndex].lStart;
		else
		{
			ReportTimeStamp("[CDSGetStart]", "ERROR:  CDS Index [" + ConvertIntToString(nIndex) + "] is Out of Range [0:" + ConvertIntToString((int)m_vCDSs.size() - 1) + "]");
		}
	}
	catch (exception ex)
	{
		cout << "ERROR [CDSGetStart] Exception Code:  " << ex.what() << "\n";
	}

	return 0;
}

////////////////////////////////////////////////////////////////////////////////
//
//  Gets the start of the CDS at the indicated index
//
////////////////////////////////////////////////////////////////////////////////
//
//  [int] nIndex :  index of the CDS to modify
//  [long] lStart:  CDS start to set
//         
////////////////////////////////////////////////////////////////////////////////

void CBase123_Catalog_Entry::CDSSetStart(int nIndex, long lStart)
{
	try
	{
		//  If the index is in range
		if ((nIndex >= 0) && (nIndex < m_vCDSs.size()))
			m_vCDSs[nIndex].lStart = lStart;
		else
		{
			ReportTimeStamp("[CDSSetStart]", "ERROR:  CDS Index [" + ConvertIntToString(nIndex) + "] is Out of Range [0:" + ConvertIntToString((int)m_vCDSs.size() - 1) + "]");
		}
	}
	catch (exception ex)
	{
		cout << "ERROR [CDSSetStart] Exception Code:  " << ex.what() << "\n";
	}
}

////////////////////////////////////////////////////////////////////////////////
//
//  Gets the stop of the CDS at the indicated index
//
////////////////////////////////////////////////////////////////////////////////
//
//  [int] nIndex:  index of the CDS to get
//  [long] lStop:  CDS stop to set
//         
////////////////////////////////////////////////////////////////////////////////

long CBase123_Catalog_Entry::CDSGetStop(int nIndex)
{
	try
	{
		//  If the index is in range
		if ((nIndex >= 0) && (nIndex < m_vCDSs.size()))
			return m_vCDSs[nIndex].lStop;
		else
		{
			ReportTimeStamp("[CDSGetStop]", "ERROR:  CDS Index [" + ConvertIntToString(nIndex) + "] is Out of Range [0:" + ConvertIntToString((int)m_vCDSs.size() - 1) + "]");
		}
	}
	catch (exception ex)
	{
		cout << "ERROR [CDSGetStop] Exception Code:  " << ex.what() << "\n";
	}

	return 0;
}

////////////////////////////////////////////////////////////////////////////////
//
//  Destructs the CBase123_Catalog_Entry class object
//
////////////////////////////////////////////////////////////////////////////////
//
//
//         
////////////////////////////////////////////////////////////////////////////////

void CBase123_Catalog_Entry::CDSSetStop(int nIndex, long lStop)
{
	try
	{
		//  If the index is in range
		if ((nIndex >= 0) && (nIndex < m_vCDSs.size()))
			m_vCDSs[nIndex].lStop = lStop;
		else
		{
			ReportTimeStamp("[CDSSetStop]", "ERROR:  CDS Index [" + ConvertIntToString(nIndex) + "] is Out of Range [0:" + ConvertIntToString((int)m_vCDSs.size() - 1) + "]");
		}
	}
	catch (exception ex)
	{
		cout << "ERROR [CDSSetStop] Exception Code:  " << ex.what() << "\n";
	}
}

////////////////////////////////////////////////////////////////////////////////
//
//  Gets the completeness of the CDS at the indicated index
//
////////////////////////////////////////////////////////////////////////////////
//
//  [int] nIndex:  index of the CDS to get
//             :  returns the completeness at the indicated index
//         
////////////////////////////////////////////////////////////////////////////////

string CBase123_Catalog_Entry::CDSGetCompleteness(int nIndex)
{
	try
	{
		//  If the index is in range
		if ((nIndex >= 0) && (nIndex < m_vCDSs.size()))
			return m_vCDSs[nIndex].strCompleteness;
		else
		{
			ReportTimeStamp("[CDSGetCompleteness]", "ERROR:  CDS Index [" + ConvertIntToString(nIndex) + "] is Out of Range [0:" + ConvertIntToString((int)m_vCDSs.size() - 1) + "]");
		}
	}
	catch (exception ex)
	{
		cout << "ERROR [CDSGetCompleteness] Exception Code:  " << ex.what() << "\n";
	}

	return "";
}

////////////////////////////////////////////////////////////////////////////////
//
//  Sets the completeness of the CDS at the indicated index
//
////////////////////////////////////////////////////////////////////////////////
//
//  [int] nIndex      :  index of the CDS to modify
//  [string] strNameID:  CDS completeness to set
//         
////////////////////////////////////////////////////////////////////////////////

void CBase123_Catalog_Entry::CDSSetCompleteness(int nIndex, string strCompleteness)
{
	try
	{
		//  If the index is in range
		if ((nIndex >= 0) && (nIndex < m_vCDSs.size()))
			m_vCDSs[nIndex].strCompleteness = strCompleteness;
		else
		{
			ReportTimeStamp("[CDSSetCompleteness]", "ERROR:  CDS Index [" + ConvertIntToString(nIndex) + "] is Out of Range [0:" + ConvertIntToString((int)m_vCDSs.size() - 1) + "]");
		}
	}
	catch (exception ex)
	{
		cout << "ERROR [CDSSetCompleteness] Exception Code:  " << ex.what() << "\n";
	}
}

////////////////////////////////////////////////////////////////////////////////
//
//  Destructs the CBase123_Catalog_Entry class object
//
////////////////////////////////////////////////////////////////////////////////
//
//
//         
////////////////////////////////////////////////////////////////////////////////

void CBase123_Catalog_Entry::ClearCDS()
{
	try
	{
		m_vCDSs.clear();
	}
	catch (exception ex)
	{
		cout << "ERROR [ClearCDS] Exception Code:  " << ex.what() << "\n";
	}
}

////////////////////////////////////////////////////////////////////////////////
//
//  Gets the entry
//
////////////////////////////////////////////////////////////////////////////////
//
//  :  returns the entry string
//         
////////////////////////////////////////////////////////////////////////////////

string CBase123_Catalog_Entry::GetEntry()
{
	//  The entry to return
	string strEntry = "";

	try
	{
		strEntry = GetEntryLine() + "|" + GetCDSLine();
	}
	catch (exception ex)
	{
		cout << "ERROR [GetEntry] Exception Code:  " << ex.what() << "\n";
	}

	return strEntry;
}

////////////////////////////////////////////////////////////////////////////////
//
//  Gets the demographics only
//
////////////////////////////////////////////////////////////////////////////////
//
//  :  returns the demographics string, if successful; else, an empty string
//         
////////////////////////////////////////////////////////////////////////////////

string CBase123_Catalog_Entry::GetDemographics()
{
	try
	{
		return GetEntryLine();
	}
	catch (exception ex)
	{
		cout << "ERROR [GetDemographics] Exception Code:  " << ex.what() << "\n";
	}

	return "";
}

////////////////////////////////////////////////////////////////////////////////
//
//  Sets the entry
//
////////////////////////////////////////////////////////////////////////////////
//
//  [string] strEntry:  entry string to parse and set
//                  :  returns true, if successful; else, false
//         
////////////////////////////////////////////////////////////////////////////////

bool CBase123_Catalog_Entry::SetEntry(string strEntry)
{
	try
	{
		//  If entry is not empty
		if (!strEntry.empty())
		{
			//  Header/CDS vector
			vector<string> vHeaderCDSs;

			//  Split header and CDS; header is required, CDS is optional
			SplitStringAllowEmptyEntries(strEntry, "|", vHeaderCDSs);

			//  If the header is properly Formatted
			if (vHeaderCDSs.size() == 2)
			{
				//  Header vector
				vector<string> vHeader;

				//  Split header
				SplitStringAllowEmptyEntries(vHeaderCDSs[0], "~", vHeader);

				//  If the header is properly Formatted
				if (vHeader.size() == 21)
				{
					m_strAccession = vHeader[0];
					m_strNameID = vHeader[1];
					ConvertLongToString(m_lLength) = vHeader[2];
					m_strMoleculeType = vHeader[3];
					m_strStrandedness = vHeader[4];
					m_strStrandednessType = vHeader[5];
					m_strStrandednessDirection = vHeader[6];
					m_strDNAIntermediate = vHeader[7];
					m_strRNAIntermediate = vHeader[8];
					m_strRetroTranscriptase = vHeader[9];
					m_strCompleteness = vHeader[10];
					m_strDate = vHeader[11];
					m_strHost = vHeader[12];
					m_strHostAge = vHeader[13];
					m_strHostGender = vHeader[14];
					m_strChromosomeSegment = vHeader[15];
					m_strSeroType = vHeader[16];
					m_strLocale = vHeader[17];
					m_strSatelliteStatus = vHeader[18];
					m_strTranscriptStatus = vHeader[19];
					m_strViralGroup = vHeader[20];

					//  If the CDS definition is not empty
					if (!vHeaderCDSs[1].empty())
					{
						//  CDSs vector
						vector<string> vCDSs;

						//  Split header
						SplitString(vHeaderCDSs[1], '~', vCDSs);

						//  If the CDS definition contains CDSs
						if (vCDSs.size() > 0)
						{
							//  Iterate through the CDSs and add
							for (int nCountCDS = 0; nCountCDS < vCDSs.size(); nCountCDS++)
							{
								//  CDS vector
								vector<string> vCDS;

								//  Split CDS
								SplitStringAllowEmptyEntries(vCDSs[nCountCDS], "^", vCDS);

								//  If the CDS definition is properly Formatted
								if (vCDS.size() == 5)
								{
									//  CDS to add
									structCDS stAdd;
									stAdd.strIsComplement = "n";
									stAdd.strNameID = "";
									stAdd.lStart = 0;
									stAdd.lStop = 0;
									stAdd.strCompleteness = "";

									//  Set CDS
									stAdd.strIsComplement = vCDS[0];
									stAdd.strNameID = vCDS[1];
									istringstream(vCDS[2]) >> stAdd.lStart;
									istringstream(vCDS[3]) >> stAdd.lStop;
									stAdd.strCompleteness = vCDS[4];

									//  Add CDS
									m_vCDSs.push_back(stAdd);
								}
							}
						}
					}
				}
				else
				{
					ReportTimeStamp("[SetEntry]", "ERROR:  Catalog Entry [" + vHeaderCDSs[0] + "] is Not Set or is Improperly Formatted - Column Count is [" + ConvertIntToString((int)vHeader.size()) + "] but Should be [21]");
				}
			}
			else
			{
				ReportTimeStamp("[SetEntry]", "ERROR:  Catalog Entry [" + strEntry.substr(0, 50) + "...] is Not Set or is Improperly Formatted - Column Count is [" + ConvertIntToString((int)vHeaderCDSs.size()) + "] but Should be [2]");
			}
		}
		else
		{
			ReportTimeStamp("[SetEntry]", "ERROR:  Catalog Entry is Empty");
		}
	}
	catch (exception ex)
	{
		cout << "ERROR [SetEntry] Exception Code:  " << ex.what() << "\n";
	}

	return false;
}

//  Implementation (private)

////////////////////////////////////////////////////////////////////////////////
//
//  Gets the [~] delmited entries string
//
////////////////////////////////////////////////////////////////////////////////
//
//  :  returns the [~] delmited entries string
//         
////////////////////////////////////////////////////////////////////////////////

string CBase123_Catalog_Entry::GetEntryLine()
{
	//  The header string to return
	string strEntries = "";

	try
	{
		strEntries = m_strAccession + "~";
		strEntries += m_strNameID + "~";
		strEntries += ConvertLongToString(m_lLength) + "~";
		strEntries += m_strMoleculeType + "~";
		strEntries += m_strStrandedness + "~";
		strEntries += m_strStrandednessType + "~";
		strEntries += m_strStrandednessDirection + "~";
		strEntries += m_strDNAIntermediate + "~";
		strEntries += m_strRNAIntermediate + "~";
		strEntries += m_strRetroTranscriptase + "~";
		strEntries += m_strCompleteness + "~";
		strEntries += m_strDate + "~";
		strEntries += m_strHost + "~";
		strEntries += m_strHostAge + "~";
		strEntries += m_strHostGender + "~";
		strEntries += m_strChromosomeSegment + "~";
		strEntries += m_strSeroType + "~";
		strEntries += m_strLocale + "~";
		strEntries += m_strSatelliteStatus + "~";
		strEntries += m_strTranscriptStatus + "~";
		strEntries += m_strViralGroup;
	}
	catch (exception ex)
	{
		cout << "ERROR [GetEntryLine] Exception Code:  [" << m_strAccession << "] " << ex.what() << "\n";
	}

	return strEntries;
}

////////////////////////////////////////////////////////////////////////////////
//
//  Gets the entry CDS starts-stops
//
////////////////////////////////////////////////////////////////////////////////
//
//  :  returns the entry CDS start-stops as a string
//         
////////////////////////////////////////////////////////////////////////////////

string CBase123_Catalog_Entry::GetCDSLine()
{
	//  The entry CDS starts-stops string to return
	string strCDSs = "";

	try
	{
		//  Cycle through CDS entries, concatenate to string
		for (int nCount = 0; nCount < m_vCDSs.size(); nCount++)
		{
			if (!strCDSs.empty())
				strCDSs += "~";

			strCDSs += m_vCDSs[nCount].strIsComplement + "^" + m_vCDSs[nCount].strNameID + "^" + ConvertLongToString(m_vCDSs[nCount].lStart) + "^" + ConvertLongToString(m_vCDSs[nCount].lStop) + "^" + m_vCDSs[nCount].strCompleteness;
		}
	}
	catch (exception ex)
	{
		cout << "ERROR [GetCDSLine] Exception Code:  " << ex.what() << "\n";
	}

	return strCDSs;
}

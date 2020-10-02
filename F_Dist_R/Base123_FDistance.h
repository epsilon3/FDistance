// Base123_FDistance.h : Performs Base123 Fofanov Distance analysis

////////////////////////////////////////////////////////////////////////////////
//
//  Performs Base123 Fofanov Distance Genomic Analysis (header); compares a foreground genome to a background genome
//      to determine genomic distance at the nucleotide polymer level; see ReportFDistanceHelp() function
//      for operational details;
//
//  Developed by Yuriy Fofanov, PhD (UTMB, yfofanov@utmb.edu)
//          and Stephen Donald Huff, PhD (Stephen.Huff.3@us.af.mil)
//          and Jared Willard (HIP Intern, Miami University)
//  Biological Informatics Group, RHDJ, 711HPW, United States Air Force Research Laboratory
//  14 October 2016
//  (All Rights Reserved)
//
////////////////////////////////////////////////////////////////////////////////

#pragma once

#include <cstdint>
#include <limits>

#include "Base123_Catalog_Entry.h"
#include "Base123_Catalog.h"

bool FilterFileForFDistance(string strInputFilePathName, bool bUseStrictFilter, int nMaxPolyLimit, string& strAcceptListFileText, string& strRejectListFileText);
bool ListFilterForFDistance(string strInputFilePathNameList, string strInputFilePathNameTransform, bool bUseStrictFilter, int nMaxPolyLimit, string strAcceptOutputListFilePathName, string strRejectOutputListFilePathName, int nMaxProcs);
bool TabulateFDistanceOutput(string& strInputPathName, string& strAccession, int nOutputCount, string& strOutputFileNameSuffix, string& strTableEntry, string& strErrorEntry);
bool ListContabulateFDistanceOutput(string strInputFilePathNameList, string strInputFilePathNameTransform, string strOutputTableFilePathName, string strCatalogFilePathName, long lMaxCatalogSize, int nOutputCount, string strOutputFileNameSuffix, string strErrorFilePathName, int nMaxProcs);
bool ListClearFDistanceOutput(string strInputFilePathNameList, string strInputFilePathNameTransform, int nOutputCount, string strOutputFileNameSuffix, string strErrorFilePathName, int nMaxProcs);
bool PerformFDistanceAnalysis(string strOutputTableFilePathName, string strBackgroundFilePathNameList, string strBackgroundInputFilePathNameTransform, string strBackgroundCatalogFilePathName, long lMaxBackgroundCatalogSize, bool bBackgroundBidirect, bool bBackgroundAllowUnknowns, string strBackgroundErrorFilePathName, string strForegroundFilePathNameList, string strForegroundInputFilePathNameTransform, string strForegroundCatalogFilePathName, long lMaxForegroundCatalogSize, bool bForegroundBidirect, bool bForegroundAllowUnknowns, string strOutputFileNameSuffix, string strForegroundErrorFilePathName, int nNMerLength, int nMaxProcs);
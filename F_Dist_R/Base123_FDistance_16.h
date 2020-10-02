// Base123_FDistance_16.h : Performs Base123 Fofanov Distance analysis

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
//  6 July 2016
//  (All Rights Reserved)
//
////////////////////////////////////////////////////////////////////////////////

#pragma once

#include <cstdint>
#include <limits>

#include "Base123_Catalog_Entry.h"
#include "Base123_Catalog.h"

uint16_t ConvertStringToSequence16(string strSequence);
string ConvertSequenceToString16(uint16_t untSequence, int nNMerLength);
bool WriteBackgroundArray16(string strOutputFilePathBaseName, int nNMerLength);
bool MutateOne16(uint16_t& untSequence);
bool MutateTwo16(uint16_t& untSequence);
bool MutateThree16(uint16_t& untSequence);
bool MutateFour16(uint16_t& untSequence);
int MutateSequence16(uint16_t& untSequence);
bool ProcessForeground16(string& strOutputFilePathName, string& strAccession, string& strSequence, bool bBidirectional, bool bForegroundAllowUnknowns, int nNMerLength, string& strOutputTableEntry);
bool MarkBackgroundSequence16(uint16_t& untSequence, int nNMerLength);
bool MutateAndMarkBackgroundSequence16(string& strSequence, int nNMerLength);
bool ProcessBackground16(string& strSequence, int nNMerLength, bool bBackgroundAllowUnknowns);
bool ProcessFDistanceList16(string strInputListFilePathName, string strInputFilePathNameTransform, CBase123_Catalog& b123Catalog, bool bBidirectional, int nNMerLength, bool bBackground, bool bAllowUnknowns, string strOutputFileNameSuffix, string strErrorFilePathName, vector<string>& vOutputTableEntries, int nMaxProcs);
bool InitializeBackground16();
bool DestroyBackground16();
bool InitializeWriteLock16();
bool DestroyWriteLock16();

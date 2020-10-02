// Base123_FDistance_32h.h : Performs Base123 Fofanov Distance analysis

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

uint32_t ConvertStringToSequence32(string strSequence);
string ConvertSequenceToString32(uint32_t untSequence, int nNMerLength);
bool WriteBackgroundArray32(string strOutputFilePathBaseName, int nNMerLength);
bool MutateOne32(uint32_t& untSequence);
bool MutateTwo32(uint32_t& untSequence);
bool MutateThree32(uint32_t& untSequence);
bool MutateFour32(uint32_t& untSequence);
bool MutateFive32(uint32_t& untSequence);
bool MutateSix32(uint32_t& untSequence);
bool MutateSeven32(uint32_t& untSequence);
bool MutateEight32(uint32_t& untSequence);
int MutateSequence32(uint32_t& untSequence);
bool ProcessForeground32(string& strOutputFilePathName, string& strAccession, string& strSequence, bool bBidirectional, bool bForegroundAllowUnknowns, int nNMerLength, string& strOutputTableEntry);
bool MarkBackgroundSequence32(uint32_t& untSequence, int nNMerLength);
bool MutateAndMarkBackgroundSequence32(string& strSequence, int nNMerLength);
bool ProcessBackground32(string& strSequence, int nNMerLength, bool bBackgroundAllowUnknowns);
bool ProcessFDistanceList32(string strInputListFilePathName, string strInputFilePathNameTransform, CBase123_Catalog& b123Catalog, bool bBidirectional, int nNMerLength, bool bBackground, bool bAllowUnknowns, string strOutputFileNameSuffix, string strErrorFilePathName, vector<string>& vOutputTableEntries, int nMaxProcs);
bool InitializeBackground32();
bool DestroyBackground32();
bool InitializeWriteLock32();
bool DestroyWriteLock32();

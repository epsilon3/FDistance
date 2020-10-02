// Base123_Utilities.h

////////////////////////////////////////////////////////////////////////////////
//
//  Base123_Utilities (header) contains general utility algorithms for inclusion in Base123 c++ code base:
//
//  Developed by Stephen Donald Huff, PhD (Stephen.Huff.3@us.af.mil)
//  Biological Informatics Group, RHDJ, 711HPW, United States Air Force Research Laboratory
//  22 June 2016
//  (All Rights Reserved)
//
////////////////////////////////////////////////////////////////////////////////

#pragma once

using namespace std;

#include <deque>

void ReportTimeStamp(string strName, string strUpdate);
bool GetFileText(string strFilePathName, stringstream& ssFileText);
bool GetFileText(string strFilePathName, string& strFileText);
bool WriteFileText(string strFilePathName, string& strFileText);
bool CreateFolderPath(string strFolderPathName);
bool ChangeWorkingFolder(string strFolderPathName);
bool IsFilePresent(string strInputFilePathName);
bool RemoveFile(string strFilePathName);
bool RenameFile(string strFromPathName, string strToPathName);
bool RenameFileByTransform(string strTransform);
bool RenameFilesByTransformSet(string strTransformSet);
bool RenameFileSetByTransformSet(vector<string> vTransformSet);
bool WriteFDistanceOutputTable(string& strOutputTableFilePathName, vector<string>& vOutputTableEntries);
string GetFDistanceOutputTableHeader();
string GetContabulatedFDistanceOutputTableHeader(int nOutputCount);
double ScoreFDistanceMutationString(string& strMutationCount, long& lTotalCount);
bool CompileFDistanceTableOutput(string& strAccession, string& strForwardOutput, string& strReverseOutput, string& strOutputTableEntry, bool bForContabulation, bool bAppendOnly);
string GetBasePath(string strInputFilePathName);
string GetBaseFileName(string strInputFilePathName);
string GetFileNameExceptLastExtension(string strInputFilePathName);
string GetFileName(string strInputFilePathName);
string GetORFOrientationFromFilePathName(string strORFFilePathName);
string TransformFilePathName(string strInputFilePathName, string strFilePathNameTransform, string strDefaultExtension);
void ParseStringToIntVector(string& strIn, char chrDelimiter, vector<int>& vTarget);
char ConvertCharacterToLowerCase(char chrToLower);
string ConvertStringToLowerCase(string strToLower);
string ConvertLongToString(long lConvert);
string ConvertDoubleToString(double dConvert);
string ConvertUnsignedInt64ToString(uint64_t nConvert);
string ConvertUnsignedInt32ToString(uint32_t nConvert);
string ConvertUnsignedInt16ToString(uint16_t nConvert);
string ConvertIntToString(int nConvert);
bool ConvertFNAtoFA(string& strFNAFileText, string& strFAFileText);
void SplitString(const string &strSplit, char chrDelimiter, vector<string> &vTarget);
void SplitString(const string &strSplit, char chrDelimiter, deque<string> &dTarget);
void SplitStringAllowEmptyEntries(const string &strSplit, string strDelimiter, vector<string> &vTarget);
bool FileTextReplace(string strInputFilePathName, string strOutputFilePathName, string strExtract, string strInsert, bool bIsCaseSensitive);
bool ListDiffVectors(vector<string>& vFrom, vector<string>& vTo, vector<string>& vDiffs);
bool VectorDiffVectors(vector<string>& vFrom, vector<string>& vTo, vector<string>& vDiffs);
bool SetDiffVectors(vector<string>& vFrom, vector<string>& vTo, vector<string>& vDiffs);
bool MultisetDiffVectors(vector<string>& vFrom, vector<string>& vTo, vector<string>& vDiffs);
string ReplaceInString(string strIn, string strExtract, string strInsert, bool bCaseSensetive);
string PadString(string strIn, string strPad, int nLength, bool bPrepend);
string ScrubDelimitedEntryString(string& strIn);
string ConvertToReverseCompliment(string strForward);
string GetAccessionFromBIGFilePathName(string strFilePathName);
string GetExtendedAccessionFromBIGFilePathName(string strFilePathName, int nMaxExtensions);
string GetAccessionFromFileHeader(string& strFileText);
string GetSequenceFromFAFile(string& strFileText);
int GetLengthFromORFFilePathName(string strORFFilePathName);
bool ParseORFFilePathName(string strORFFilePathName, string& strOrientation, long& lStart, long& lStop, long& lLength);
string ConvertORFToFAFilePathName(string strORFFilePathName);
string GetErrorMessage(int nError);

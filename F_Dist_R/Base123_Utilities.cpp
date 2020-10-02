// Base123_Utilities.cpp

////////////////////////////////////////////////////////////////////////////////
//
//  Contains general utility algorithms for inclusion in Base123 c++ code base:
//
//  Developed by Stephen Donald Huff, PhD (Stephen.Huff.3@us.af.mil)
//  Biological Informatics Group, RHDJ, 711HPW, United States Air Force Research Laboratory
//  22 June 2016
//  (All Rights Reserved)
//
////////////////////////////////////////////////////////////////////////////////

#include "F_Dist_R.h"
#include "Base123_Utilities.h"

#include <iostream>
#include <fstream>
#include <sstream>
#include <ctime>
#include <locale>
#include <sys/types.h>
#include <sys/stat.h>
#include <unordered_set>
#include <omp.h>
#include <deque>
#include <iterator>
#include <set>
#include <algorithm>

#ifdef _WIN64
	#include <direct.h>
#else
	#include <unistd.h>
	#include <limits.h>
#endif

////////////////////////////////////////////////////////////////////////////////
//
//  Reports a time-stamped update to the console; format is:
//      Name: Update: Year-Month-Day_of_Month Hour:Minute:Second
//
////////////////////////////////////////////////////////////////////////////////
//
//  [string] strName  :  name of the update object
//  [string] strUpdate:  update descripton
//         
////////////////////////////////////////////////////////////////////////////////

void ReportTimeStamp(string strName, string strUpdate)
{
	try
	{
		#ifdef _WIN64
			//  Current time
			const time_t tStart = time(0);
			//  Current time structure
			struct tm tmStart;

			localtime_s(&tmStart, &tStart);

			//  Report
			cout << strName << ": " << strUpdate << ": " << (1900 + tmStart.tm_year) << "-" << tmStart.tm_mon + 1 << "-" << tmStart.tm_mday << " " << tmStart.tm_hour << ":" << tmStart.tm_min << ":" << tmStart.tm_sec << "\n";
		#else
			//  Current time
			const time_t tStart = time(0);
			//  Current time structure
			struct tm* tmStart;

			tmStart = localtime(&tStart);

			//  Report
			cout << strName << ": " << strUpdate << ": " << (1900 + tmStart->tm_year) << "-" << tmStart->tm_mon + 1 << "-" << tmStart->tm_mday << " " << tmStart->tm_hour << ":" << tmStart->tm_min << ":" << tmStart->tm_sec << "\n";
		#endif

	}
	catch (exception ex)
	{
		cout << "ERROR [ReportTimeStamp] Exception Code:  " << ex.what() << "\n";
	}

	return;
}

////////////////////////////////////////////////////////////////////////////////
//
//  Gets a string representing the contents of a text file;
//         
////////////////////////////////////////////////////////////////////////////////
//
//  [string] strFilePathName  :  file path name of the file to read and return;
//                                      NOTE:  strips carriage-return characters
//  [istringstream&] isileText:  file text to fill
//                            :  returns true, if successful; else, false
//         
////////////////////////////////////////////////////////////////////////////////

bool GetFileText(string strFilePathName, stringstream& ssFileText)
{
	//  The input file
	ifstream iFile;

	try
	{
		//  If the file path name is not empty, then open the file
		if (!strFilePathName.empty())
		{
			//  Open the file
			iFile.open(strFilePathName.c_str(), ios::in);

			//  If the file is open, then read the file
			if (iFile.is_open())
			{
				//  Read file
				ssFileText << iFile.rdbuf();

				//  Close the file
				iFile.close();

				return true;
			}
			else
			{
				ReportTimeStamp("[GetFileText]", "ERROR:  Input File [" + strFilePathName + "] Open Failed");
			}
		}
		else
		{
			ReportTimeStamp("[GetFileText]", "ERROR:  Input File Path Name is Empty");
		}
	}
	catch (exception ex)
	{
		cout << "ERROR [GetFileText] Exception Code:  " << ex.what() << "\n";
	}

	return false;
}

////////////////////////////////////////////////////////////////////////////////
//
//  Gets a string representing the contents of a text file;
//         
////////////////////////////////////////////////////////////////////////////////
//
//  [string] strFilePathName:  file path name of the file to read and return;
//                                      NOTE:  strips carriage-return characters
//  [string] strFileText    :  file text to fill
//                         :  returns true, if successful; else, false
//         
////////////////////////////////////////////////////////////////////////////////

bool GetFileText(string strFilePathName, string& strFileText)
{
	//  The input file
	ifstream iFile;
	//  The read character
	char chrGet;

	try
	{
		//  Reset file text
		strFileText = "";

		//  If the file path name is not empty, then open the file
		if (!strFilePathName.empty())
		{
			//  Open the file
			iFile.open(strFilePathName.c_str(), ios::in);

			//  If the file is open, then read the file
			if (iFile.is_open())
			{
				//  Read the file
				while (iFile.get(chrGet))
				{
					//  Strip carriage-return chracters, else concatenate return string
					if (chrGet != '\r')
						strFileText += chrGet;
				}

				//  Close the file
				iFile.close();

				return true;
			}
			else
			{
				ReportTimeStamp("[GetFileText]", "ERROR:  Input File [" + strFilePathName + "] Open Failed");
			}
		}
		else
		{
			ReportTimeStamp("[GetFileText]", "ERROR:  Input File Path Name is Empty");
		}
	}
	catch (exception ex)
	{
		cout << "ERROR [GetFileText] Exception Code:  " << ex.what() << "\n";
	}

	return false;
}

////////////////////////////////////////////////////////////////////////////////
//
//  Writes a string to a text file;
//         
////////////////////////////////////////////////////////////////////////////////
//
//  [string] strFilePathName:  file path name of the file to write;
//  [string] strFileText    :  text to write;
//                         :  returns true, if successful; else, false
//         
////////////////////////////////////////////////////////////////////////////////

bool WriteFileText(string strFilePathName, string& strFileText)
{
	//  The output file
	ofstream oFile;

	try
	{
		//  If the file path name is not empty, then open the file
		if (!strFilePathName.empty())
		{
			//  Open the file
			oFile.open(strFilePathName.c_str(), ios::out);

			//  If the file is open, then write to the file
			if (oFile.is_open())
			{
				//  Write to the file
				//oFile << strFileText;
				oFile.write(strFileText.c_str(), sizeof(char)*strFileText.length());

				#ifdef _WIN64
				#else

					//  Set the umask
					umask(0007);

					chmod(strFilePathName.c_str(), S_IRUSR | S_IWUSR | S_IRGRP | S_IWGRP);

				#endif

				//  Close the file
				oFile.close();

				return true;
			}
			else
			{
				ReportTimeStamp("[WriteFileText]", "ERROR:  Output File [" + strFilePathName + "] Create/Open Failed");
			}
		}
		else
		{
			ReportTimeStamp("[WriteFileText]", "ERROR:  Output File Path Name is Empty");
		}
	}
	catch (exception ex)
	{
		cout << "ERROR [WriteFileText] Exception Code:  " << ex.what() << "\n";
	}

	return false;
}

////////////////////////////////////////////////////////////////////////////////
//
//  Creates the requested folder path
//
////////////////////////////////////////////////////////////////////////////////
//
//  [string] strFolderPathName:  folder path name to be created
//         
////////////////////////////////////////////////////////////////////////////////

bool CreateFolderPath(string strFolderPathName)
{
	try
	{
		//  If folder path name is not empty, then create the folder path
		if (!strFolderPathName.empty())
		{
			struct stat stGet = { 0 };

			if (stat(strFolderPathName.c_str(), &stGet) == -1)
			{
				#ifdef _WIN64
					_mkdir(strFolderPathName.c_str());
				#else
					//  Set the umask
					umask(0007);

					//  S_IRWXU = Read-Write-Execute-Search (Owner)
					//  S_IRWXG = Read-Write-Execute-Search (Group)
					mkdir(strFolderPathName.c_str(), S_IRWXU | S_IRWXG);
				#endif
			}

			return true;
		}
		else
		{
			ReportTimeStamp("[CreateFolderPath]", "ERROR:  Folder Path Name is Empty");
		}
	}
	catch (exception ex)
	{
		cout << "ERROR [CreateFolderPath] Exception Code:  " << ex.what() << "\n";
	}

	return false;
}

////////////////////////////////////////////////////////////////////////////////
//
//  Changes working focus to the requested folder path
//
////////////////////////////////////////////////////////////////////////////////
//
//  [string] strFolderPathName:  folder path name to be created
//         
////////////////////////////////////////////////////////////////////////////////

bool ChangeWorkingFolder(string strFolderPathName)
{
	try
	{
		//  If folder path name is not empty, then create the folder path
		if (!strFolderPathName.empty())
		{
			#ifdef _WIN64
			#else
				if (chdir(strFolderPathName.c_str()) == -1)
					ReportTimeStamp("[ChangeWorkingFolder]", "ERROR:  Working Folder Change Failed [" + GetErrorMessage(errno) + "]");
				else
				{
					/*char* pCWD = NULL;
					
					pCWD = new char[PATH_MAX];
					if (pCWD != NULL)
					{
						pCWD = getcwd(pCWD, PATH_MAX);
						string strCWD(pCWD);
						ReportTimeStamp("[ChangeWorkingFolder]", "NOTE:  Working Folder is [" + strCWD + "]");
					}*/

					return true;
				}
			#endif
		}
		else
		{
			ReportTimeStamp("[ChangeWorkingFolder]", "ERROR:  Folder Path Name is Empty");
		}
	}
	catch (exception ex)
	{
		cout << "ERROR [ChangeWorkingFolder] Exception Code:  " << ex.what() << "\n";
	}

	return false;
}

////////////////////////////////////////////////////////////////////////////////
//
//  Checks for presence of file
//         
////////////////////////////////////////////////////////////////////////////////
//
//  [string] strInputFilePathName:  file path name to check
//                               :  returns true, if file is present; else, false
//         
////////////////////////////////////////////////////////////////////////////////

bool IsFilePresent(string strInputFilePathName)
{
	try
	{
		//  If input file path name is not empty
		if (!strInputFilePathName.empty())
		{
			#ifdef _WIN64
			#else
				//  If file exists
				if (access(strInputFilePathName.c_str(), 0) == 0)
					return true;
			#endif
		}
	}
	catch (exception ex)
	{
		cout << "ERROR [IsFilePresent] Exception Code:  " << ex.what() << "\n";
	}

	return false;
}

////////////////////////////////////////////////////////////////////////////////
//
//  Removes file
//         
////////////////////////////////////////////////////////////////////////////////
//
//  [string] strFilePathName:  file path name to remove
//                          :  returns true, if successful; else, false
//         
////////////////////////////////////////////////////////////////////////////////

bool RemoveFile(string strFilePathName)
{
	try
	{
		//  If file path name is not empty
		if (!strFilePathName.empty())
		{
			#ifdef _WIN64
			#else
				if (IsFilePresent(strFilePathName))
				{
					if (remove(strFilePathName.c_str()) == 0)
						return true;
				}
			#endif
		}
	}
	catch (exception ex)
	{
		cout << "ERROR [RemoveFile] Exception Code:  " << ex.what() << "\n";
	}

	return false;
}

////////////////////////////////////////////////////////////////////////////////
//
//  Renames files
//         
////////////////////////////////////////////////////////////////////////////////
//
//  [string] strFromPathName:  target file path name to rename
//  [string] strToPathName  :  destination file path name
//                         :  returns true, if successful; else, false
//         
////////////////////////////////////////////////////////////////////////////////

bool RenameFile(string strFromPathName, string strToPathName)
{
	try
	{
		//  If from and to path name are not empty
		if ((!strFromPathName.empty()) && (!strToPathName.empty()))
		{
			//ReportTimeStamp("[RenameFile]", "NOTE:  FROM [" + strFromPathName + "] to [" + strToPathName + "]");

			#ifdef _WIN64
			#else
				//  If file exists
				if (IsFilePresent(strFromPathName))
				{
					//  Rename
					if (rename(strFromPathName.c_str(), strToPathName.c_str()) == 0)
						return true;
					else
					{
						ReportTimeStamp("[RenameFile]", "ERROR:  Rename [" + strFromPathName + "] to [" + strToPathName + "] Failed [" + GetErrorMessage(errno) + "]");

						return false;
					}
				}
			#endif
		}
	}
	catch (exception ex)
	{
		cout << "ERROR [RenameFile] Exception Code:  " << ex.what() << "\n";
	}

	return false;
}

////////////////////////////////////////////////////////////////////////////////
//
//  Rename file by transform
//         
////////////////////////////////////////////////////////////////////////////////
//
//  [string] strTransform:  target to destination file path name to rename; entries delimited in~out
//                      :  returns true, if successful; else, false
//         
////////////////////////////////////////////////////////////////////////////////

bool RenameFileByTransform(string strTransform)
{
	//  File path names
	vector<string> vFilePathNames;

	try
	{
		//  If transform is not empty
		if (!strTransform.empty())
		{
			//  Split fromm and to file path names
			SplitString(strTransform, '~', vFilePathNames);

			//  If from and to file path name are set
			if (vFilePathNames.size() == 2)
			{
				//cout << "[RenameFileByTransform] NOTE:  FAILED " << vFilePathNames[0] << " >> " << vFilePathNames[1] << "\n";

				//  Rename file
				return RenameFile(vFilePathNames[0], vFilePathNames[1]);
			}
		}
	}
	catch (exception ex)
	{
		cout << "ERROR [RenameFileByTransform] Exception Code:  " << ex.what() << "\n";
	}

	return false;
}

////////////////////////////////////////////////////////////////////////////////
//
//  Renames files by transform set
//         
////////////////////////////////////////////////////////////////////////////////
//
//  [string] strTransformSet:  target to destination file path names to rename; entries delimited in~out^in~out^...
//                         :  returns true, if successful; else, false
//         
////////////////////////////////////////////////////////////////////////////////

bool RenameFilesByTransformSet(string strTransformSet)
{
	//  Return state
	bool bReturnState = true;
	//  Transforms
	vector<string> vTransforms;

	try
	{
		//  If transforms is not empty
		if (!strTransformSet.empty())
		{
			//  Split transforms
			SplitString(strTransformSet, '&', vTransforms);

			//  If transforms container is not empty
			if (vTransforms.size() > 0)
			{
				//  Iterate transforms, rename files
				for (int nCount = 0; nCount < (int)vTransforms.size(); nCount++)
				{
					if (!RenameFileByTransform(vTransforms[nCount]))
					{
						//cout << "[RenameFilesByTransformSet] NOTE:  FAILED " << vTransforms[nCount] << "\n";

						bReturnState = false;
					}
					else
					{
						//cout << "[RenameFilesByTransformSet] NOTE:  SUCCEEDED " << vTransforms[nCount] << "\n";

					}
				}
			}
		}
	}
	catch (exception ex)
	{
		cout << "ERROR [RenameFilesByTransformSet] Exception Code:  " << ex.what() << "\n";
	}

	return bReturnState;
}

////////////////////////////////////////////////////////////////////////////////
//
//  Renames file set by transform set
//         
////////////////////////////////////////////////////////////////////////////////
//
//  [vector<string>] vTransformSet:  target to destination file path names to rename; vector entries delimited in~out^in~out^...
//                               :  returns true, if successful; else, false
//         
////////////////////////////////////////////////////////////////////////////////

bool RenameFileSetByTransformSet(vector<string> vTransformSet)
{
	//  Return state
	bool bReturnState = true;

	try
	{
		//  If transforms set is not empty
		if (vTransformSet.size() > 0)
		{
			//  Iterate transforms set
			for (int nCount = 0; nCount < (int)vTransformSet.size(); nCount++)
			{
				if (!RenameFilesByTransformSet(vTransformSet[nCount]))
					bReturnState = false;
			}
		}
	}
	catch (exception ex)
	{
		cout << "ERROR [RenameFileSetByTransformSet] Exception Code:  " << ex.what() << "\n";
	}

	return bReturnState;
}

////////////////////////////////////////////////////////////////////////////////
//
//  Writes the F-Distance output table
//         
////////////////////////////////////////////////////////////////////////////////
//
//  [string&] strOutputTableFilePathName  :  output table file path name
//  [vector<string>&] vOutputTableEntries:  output table file text to write
//                                       :  returns true, if successful; else, false
//         
////////////////////////////////////////////////////////////////////////////////

bool WriteFDistanceOutputTable(string& strOutputTableFilePathName, vector<string>& vOutputTableEntries)
{
	//  File text
	string strFileText = "";

	try
	{
		//  If output table file path name is not empty
		if (!strOutputTableFilePathName.empty())
		{
			//  If output table file path name is not empty
			if (vOutputTableEntries.size() > 0)
			{
				//  Add header
				strFileText = GetFDistanceOutputTableHeader();

				//  Iterate entries and concatenate file text
				for (long lCount = 0; lCount < vOutputTableEntries.size(); lCount++)
				{
					if (!vOutputTableEntries[lCount].empty())
						strFileText += vOutputTableEntries[lCount] + "\n";
				}

				//  Write output table
				return WriteFileText(strOutputTableFilePathName, strFileText);
			}
			else
			{
				ReportTimeStamp("[WriteFDistanceOutputTable]", "ERROR:  F-Distance Output Table File Text is Empty");
			}
		}
		else
		{
			ReportTimeStamp("[WriteFDistanceOutputTable]", "ERROR:  F-Distance Output Table File Path Name is Empty");
		}
	}
	catch (exception ex)
	{
		cout << "ERROR [WriteFDistanceOutputTable] Exception Code:  " << ex.what() << "\n";
	}

	return false;
}

////////////////////////////////////////////////////////////////////////////////
//
//  Gets the F-Distance output table header
//         
////////////////////////////////////////////////////////////////////////////////
//
//  :  returns output table header, if successful; else, empty string
//         
////////////////////////////////////////////////////////////////////////////////

string GetFDistanceOutputTableHeader()
{
	//  Header to return
	string strHeader = "";

	try
	{
		strHeader += "Accession~Length~Forward Mutation Total~Forward F-Distance Score~Reverse Mutation Total~Reverse F-Distance Score\n";

		return strHeader;
	}
	catch (exception ex)
	{
		cout << "ERROR [GetFDistanceOutputTableHeader] Exception Code:  " << ex.what() << "\n";
	}

	return strHeader;
}

////////////////////////////////////////////////////////////////////////////////
//
//  Gets the contabulated F-Distance output table header
//         
////////////////////////////////////////////////////////////////////////////////
//
//  [int] nOutputCount:  original Shuffler output file count
//                   :  returns output table header, if successful; else, empty string
//         
////////////////////////////////////////////////////////////////////////////////

string GetContabulatedFDistanceOutputTableHeader(int nOutputCount)
{
	//  Header to return
	string strHeader = "";

	try
	{
		strHeader += "~Length (F-Distance)~Forward Mutation Total (gen)~Forward F-Distance Score (gen)~Reverse Mutation Total (gen)~Reverse F-Distance Score (gen)";

		//  Append header entries for shuffle file output
		for (int nCount = 0; nCount < nOutputCount; nCount++)
			strHeader += "~Forward Mutation Total (sh_" + ConvertIntToString(nCount) + ")~Forward F-Distance Score (sh_" + ConvertIntToString(nCount) + ")~Reverse Mutation Total (sh_" + ConvertIntToString(nCount) + ")~Reverse F-Distance Score (sh_" + ConvertIntToString(nCount) + ")";

		strHeader += "\n";

		return strHeader;
	}
	catch (exception ex)
	{
		cout << "ERROR [GetContabulatedFDistanceOutputTableHeader] Exception Code:  " << ex.what() << "\n";
	}

	return strHeader;
}

////////////////////////////////////////////////////////////////////////////////
//
//  Scores an F-Distance mutation count string
//         
////////////////////////////////////////////////////////////////////////////////
//
//  [string&] strMutationCount:  F-Distance mutation count string
//  [long&] lTotalCount       :  total mutation count to return
//                           :  returns F-Distance mutation count score [double], if successful; else, -1
//         
////////////////////////////////////////////////////////////////////////////////

double ScoreFDistanceMutationString(string& strMutationCount, long& lTotalCount)
{
	//  Score to return
	double dScore = -1;
	//  Get score
	int nGetScore = 0;

	try
	{
		//  Zero total count
		lTotalCount = 0;

		//  If mutation count string is not empty
		if (!strMutationCount.empty())
		{
			//  Iterate and total
			for (long lCount = 0; (lCount < strMutationCount.length()) && (lCount < strMutationCount.length()); lCount++)
			{
				//  Get the score at this position
				stringstream(strMutationCount.substr(lCount, 1)) >> nGetScore;
				//  Add to total
				lTotalCount += nGetScore;
			}

			//  Score
			dScore = (double)((double)lTotalCount / (double)strMutationCount.length());
		}
		else
		{
			ReportTimeStamp("[ScoreFDistanceMutationString]", "ERROR:  Mutation Count Container is Empty");
		}
	}
	catch (exception ex)
	{
		cout << "ERROR [ScoreFDistanceMutationString] Exception Code:  " << ex.what() << "\n";
	}

	return dScore;
}

////////////////////////////////////////////////////////////////////////////////
//
//  Gets the F-Distance output table header
//         
////////////////////////////////////////////////////////////////////////////////
//
//  [string&] strAccession       :  sequence accession
//  [string&] strForwardOutput   :  forawrd F-Distance mutation
//  [string&] strReverseOutput   :  reverse F-Distance mutation scores
//  [string&] strOutputTableEntry:  output table file text
//  [bool] bForContabulation     :  if true, compilation is for contabulation output; else, compilation is for stand-alone output
//  [bool] bAppendOnly           :  if true, append only the score to the table entry text; else, set entire table entry text
//                              :  returns output table header, if successful; else, empty string
//         
////////////////////////////////////////////////////////////////////////////////

bool CompileFDistanceTableOutput(string& strAccession, string& strForwardOutput, string& strReverseOutput, string& strOutputTableEntry, bool bForContabulation, bool bAppendOnly)
{
	//  Forward total
	long lTotalForward = 0;
	//  Reverse total
	long lTotalReverse = 0;
	//  Forward score
	double dScoreForward = 0;
	//  Reverse score
	double dScoreReverse = 0;
	//  Get score
	int nGetScore = 0;

	try
	{
		//  If accession is not empty
		if (!strAccession.empty())
		{
			//  Lengths must be equal
			if (((!strForwardOutput.empty()) && (!strReverseOutput.empty()) && (strForwardOutput.length() == strReverseOutput.length())) || ((!strForwardOutput.empty()) || (!strReverseOutput.empty())))
			{
				//  Score
				if (!strForwardOutput.empty())
					dScoreForward = ScoreFDistanceMutationString(strForwardOutput, lTotalForward);
				if(!strReverseOutput.empty())
					dScoreReverse = ScoreFDistanceMutationString(strReverseOutput, lTotalReverse);

				//  Append scores to existing score entry
				if (bAppendOnly)
				{
					strOutputTableEntry += "~" + ConvertLongToString(lTotalForward) + "~";
					strOutputTableEntry += ConvertDoubleToString(dScoreForward) + "~";
					strOutputTableEntry += ConvertLongToString(lTotalReverse) + "~";
					strOutputTableEntry += ConvertDoubleToString(dScoreReverse);
				}
				//  Create entire score entry
				else
				{
					//  Contabulation includes accession from catalog, append
					if(bForContabulation)
						strOutputTableEntry += "~";
					//  Standalone requires accession, set new entry
					else
						strOutputTableEntry = strAccession + "~";

					//  Bot require F-Distance length to double-check NCBI length
					if (!strForwardOutput.empty())
						strOutputTableEntry += ConvertLongToString((long)strForwardOutput.length()) + "~";
					else if (!strReverseOutput.empty())
						strOutputTableEntry += ConvertLongToString((long)strReverseOutput.length()) + "~";

					strOutputTableEntry += ConvertLongToString(lTotalForward) + "~";
					strOutputTableEntry += ConvertDoubleToString(dScoreForward) + "~";
					strOutputTableEntry += ConvertLongToString(lTotalReverse) + "~";
					strOutputTableEntry += ConvertDoubleToString(dScoreReverse);					
				}

				return true;
			}
			else
			{
				ReportTimeStamp("[CompileFDistanceTableOutput]", "ERROR:  Forward (and/or Reverse) Counts are Not Properly Formatted");
			}
		}
		else
		{
			ReportTimeStamp("[CompileFDistanceTableOutput]", "ERROR:  Accession is Empty");
		}
	}
	catch (exception ex)
	{
		cout << "ERROR [CompileFDistanceTableOutput] Exception Code:  " << ex.what() << "\n";
	}

	return false;
}

////////////////////////////////////////////////////////////////////////////////
//
//  Parses input file path name for base path
//         
////////////////////////////////////////////////////////////////////////////////
//
//  [string] strInputFilePathName:  file path to parse
//                              :  returns the base path, if successful; else, an empty string
//         
////////////////////////////////////////////////////////////////////////////////

string GetBasePath(string strInputFilePathName)
{
	//  Base path
	string strBasePath = "";
	//  Path delimiter character
	char chrPathDelimiter = strPathDelimiter;
	//  File path name parts
	vector<string> vFilePathNameParts;

	try
	{
		//  If input file path name is not empty
		if (!strInputFilePathName.empty())
		{
			//  Split file path name parts
			SplitString(strInputFilePathName, chrPathDelimiter, vFilePathNameParts);

			//  If file path name parts is not empty
			if (vFilePathNameParts.size() >= 1)
			{
				//  Iterate file path name parts
				for (int nCountPathParts = 0; nCountPathParts < (int)vFilePathNameParts.size() - 1; nCountPathParts++)
					strBasePath += chrPathDelimiter + vFilePathNameParts[nCountPathParts];
			}
		}
	}
	catch (exception ex)
	{
		cout << "ERROR [GetBasePath] Exception Code:  " << ex.what() << "\n";
	}

	return strBasePath;
}

////////////////////////////////////////////////////////////////////////////////
//
//  Parses for file name, only, no extensions
//         
////////////////////////////////////////////////////////////////////////////////
//
//  [string] strInputFilePathName:  file path name to parse
//                              :  returns the file name without extension, if successful; else, an empty string
//         
////////////////////////////////////////////////////////////////////////////////

string GetBaseFileName(string strInputFilePathName)
{
	//  File name to return
	string strFileName = "";
	//  Path delimiter character
	char chrPathDelimiter = strPathDelimiter;
	//  File path name parts
	vector<string> vFilePathNameParts;
	//  File name parts
	vector<string> vFileNameParts;

	try
	{
		//  If input file path name is not empty
		if (!strInputFilePathName.empty())
		{
			//  Split file path name parts
			SplitString(strInputFilePathName, chrPathDelimiter, vFilePathNameParts);

			//  If file path name parts is not empty
			if (vFilePathNameParts.size() >= 1)
			{
				//  Split file name parts
				SplitString(vFilePathNameParts[vFilePathNameParts.size() - 1], '.', vFileNameParts);

				//  If file name parts is not empty
				if (vFileNameParts.size() > 0)
					strFileName = vFileNameParts[0];
			}
		}
	}
	catch (exception ex)
	{
		cout << "ERROR [GetBaseFileName] Exception Code:  " << ex.what() << "\n";
	}

	return strFileName;
}

////////////////////////////////////////////////////////////////////////////////
//
//  Parses for input file name without final extension
//         
////////////////////////////////////////////////////////////////////////////////
//
//  [string] strInputFilePathName:  file path name to parse
//                               :  returns the file name without extension, if successful; else, an empty string
//         
////////////////////////////////////////////////////////////////////////////////

string GetFileNameExceptLastExtension(string strInputFilePathName)
{
	//  File name to return
	string strFileName = "";
	//  Path delimiter character
	char chrPathDelimiter = strPathDelimiter;
	//  File path name parts
	vector<string> vFilePathNameParts;
	//  File name parts
	vector<string> vFileNameParts;

	try
	{
		//  If input file path name is not empty
		if (!strInputFilePathName.empty())
		{
			//  Split file path name parts
			SplitString(strInputFilePathName, chrPathDelimiter, vFilePathNameParts);

			//  If file path name parts is not empty
			if (vFilePathNameParts.size() >= 1)
			{
				//  Split file name parts
				SplitString(vFilePathNameParts[vFilePathNameParts.size() - 1], '.', vFileNameParts);

				//  If file name parts is not empty
				if (vFileNameParts.size() > 0)
				{
					//  Iterate parts, concatenate all save final part
					for (int nCount = 0; nCount < vFileNameParts.size() - 1; nCount++)
					{
						//  If file name is not empty, prepend '.'
						if (!strFileName.empty())
							strFileName += ".";
						strFileName += vFileNameParts[nCount];
					}
				}
			}
		}
	}
	catch (exception ex)
	{
		cout << "ERROR [GetFileNameExceptLastExtension] Exception Code:  " << ex.what() << "\n";
	}

	return strFileName;
}

////////////////////////////////////////////////////////////////////////////////
//
//  Parses input file path name for file name with extension
//         
////////////////////////////////////////////////////////////////////////////////
//
//  [string] strInputFilePathName:  file path name to parse
//                              :  returns the file name without extension, if successful; else, an empty string
//         
////////////////////////////////////////////////////////////////////////////////

string GetFileName(string strInputFilePathName)
{
	//  File name to return
	string strFileName = "";
	//  Path delimiter character
	char chrPathDelimiter = strPathDelimiter;
	//  File path name parts
	vector<string> vFilePathNameParts;
	//  File name parts
	vector<string> vFileNameParts;

	try
	{
		//  If input file path name is not empty
		if (!strInputFilePathName.empty())
		{
			//  Split file path name parts
			SplitString(strInputFilePathName, chrPathDelimiter, vFilePathNameParts);

			//  If file path name parts is not empty
			if (vFilePathNameParts.size() >= 1)
			{
				//  Split file name parts
				SplitString(vFilePathNameParts[vFilePathNameParts.size() - 1], '.', vFileNameParts);

				//  If file name parts is not empty
				if (vFileNameParts.size() > 0)
				{
					//  Iterate parts, concatenate all
					for (int nCount = 0; nCount < vFileNameParts.size(); nCount++)
					{
						//  If file name is not empty, prepend '.'
						if (!strFileName.empty())
							strFileName += ".";
						strFileName += vFileNameParts[nCount];
					}
				}
			}
		}
	}
	catch (exception ex)
	{
		cout << "ERROR [GetFileName] Exception Code:  " << ex.what() << "\n";
	}

	return strFileName;
}

////////////////////////////////////////////////////////////////////////////////
//
//  Gets ORF orientation (f/r) from .ORF file path name
//         
////////////////////////////////////////////////////////////////////////////////
//
//  [string] strORFFilePathName:  .ORF file path name to parse
//                             :  returns ORF orientation, if successful; else, an empty string
//         
////////////////////////////////////////////////////////////////////////////////

string GetORFOrientationFromFilePathName(string strORFFilePathName)
{
	//  ORF Orientation to return
	string strORFOrientation = "";

	try
	{
		//  If file name is not empty
		if (!strORFFilePathName.empty())
		{
			if ((ConvertStringToLowerCase(strORFFilePathName).find(".f1_") != string::npos) || (ConvertStringToLowerCase(strORFFilePathName).find(".f2_") != string::npos) ||
				(ConvertStringToLowerCase(strORFFilePathName).find(".f3_") != string::npos))
				strORFOrientation = "f";
			else if ((ConvertStringToLowerCase(strORFFilePathName).find(".r1_") != string::npos) || (ConvertStringToLowerCase(strORFFilePathName).find(".r2_") != string::npos) ||
				(ConvertStringToLowerCase(strORFFilePathName).find(".r3_") != string::npos))
				strORFOrientation = "r";
		}
		else
		{
			ReportTimeStamp("[GetORFOrientationFromFilePathName]", "ERROR:  .ORF File Path Name is Empty");
		}
	}
	catch (exception ex)
	{
		cout << "ERROR [GetORFOrientationFromFilePathName] Exception Code:  " << ex.what() << "\n";
	}

	return strORFOrientation;
}

////////////////////////////////////////////////////////////////////////////////
//
//  Transforms a file name into the indicated file path name (getting BIG accession number from file name and replacing into transform)
//         
////////////////////////////////////////////////////////////////////////////////
//
//  [string] strInputFilePathName    :  string to parse
//  [string] strFilePathNameTransform:  element delimiter character
//  [string] strDefaultExtension     :  default file extension
//                                   :  returns tranformed file path name, if successful; else, an empty string
//         
////////////////////////////////////////////////////////////////////////////////

string TransformFilePathName(string strInputFilePathName, string strFilePathNameTransform, string strDefaultExtension)
{
	//  Transformed file path name to return
	string strFilePathName = "";
	//  Accession number
	string strAccession = "";
	//  Path delimiter
	char chrPathDelimiter = strPathDelimiter;

	try
	{
		//  If file name is not empty
		if (!strInputFilePathName.empty())
		{
			//  If transform string is not empty
			if (!strFilePathNameTransform.empty())
			{
				//  Get BIG accession
				strAccession = GetAccessionFromBIGFilePathName(strInputFilePathName);

				//  If accession is not empty
				if (!strAccession.empty())
				{
					strFilePathName = ReplaceInString(strFilePathNameTransform, "^BIG_ACCESSION^", strAccession, true);
					strFilePathName = ReplaceInString(strFilePathName, "^BASE_FILE_PATH_NAME^", GetBasePath(strInputFilePathName) + chrPathDelimiter + GetFileNameExceptLastExtension(strInputFilePathName), true);
					strFilePathName = ReplaceInString(strFilePathName, "^BASE_PATH^", GetBasePath(strInputFilePathName), true);
					strFilePathName = ReplaceInString(strFilePathName, "^BASE_FILE_NAME^", GetBaseFileName(strInputFilePathName), true);
					strFilePathName = ReplaceInString(strFilePathName, "^FILE_NAME^", GetFileName(strInputFilePathName), true);
					strFilePathName = ReplaceInString(strFilePathName, "^FNELE^", GetFileNameExceptLastExtension(strInputFilePathName), true);
					strFilePathName = ReplaceInString(strFilePathName, "^DEFAULT_EXTENSION^", strDefaultExtension, true);
				}
				else
				{
					ReportTimeStamp("[TransformFilePathName]", "ERROR:  Accession is Empty");
				}
			}
			else
			{
				ReportTimeStamp("[TransformFilePathName]", "ERROR:  File Path Name Transform is Empty");
			}
		}
		else
		{
			ReportTimeStamp("[TransformFilePathName]", "ERROR:  File Name is Empty");
		}
	}
	catch (exception ex)
	{
		cout << "ERROR [TransformFilePathName] Exception Code:  " << ex.what() << "\n";
	}

	return strFilePathName;
}

////////////////////////////////////////////////////////////////////////////////
//
//  Parses a delimited string to a vector of integers
//         
////////////////////////////////////////////////////////////////////////////////
//
//  [string&] strIn    :  string to parse
//  [char] chrDelimiter:  element delimiter character
//  [vector<int>] vIn  :  vector<int> to fill
//         
////////////////////////////////////////////////////////////////////////////////

void ParseStringToIntVector(string& strIn, char chrDelimiter, vector<int>& vTarget)
{
	//  The temporary vector<string>
	vector<string> vSplit;

	try
	{
		//  Clear input vector
		vTarget.clear();

		//  If input string is not empty
		if (!strIn.empty())
		{
			//  If delmiter is not empty
			if (chrDelimiter != '\0')
			{
				//  Split the string to vector<string>
				SplitString(strIn, chrDelimiter, vSplit);

				//  If vector is set
				if (vSplit.size() > 0)
				{
					//  Push-back ints into vTarget
					for (int nCount = 0; nCount < vSplit.size(); nCount++)
					{
						//  Integer conversion target
						int nGet = 0;

						//  Set maximum processor count
						istringstream(vSplit[nCount]) >> nGet;

						//  Push-back the integer
						vTarget.push_back(nGet);
					}
				}
				//else
				//{
				//	ReportTimeStamp("[ParseStringToIntVector]", "ERROR:  Input is Not Set");
				//}
			}
			else
			{
				ReportTimeStamp("[ParseStringToIntVector]", "ERROR:  Delimiter is Empty");
			}
		}
		//else
		//{
		//	ReportTimeStamp("[ParseStringToIntVector]", "ERROR:  Input is Empty");
		//}
	}
	catch (exception ex)
	{
		cout << "ERROR [ParseStringToIntVector] Exception Code:  " << ex.what() << "\n";
	}

	return;
}

////////////////////////////////////////////////////////////////////////////////
//
//  Converts a character to lower-case in a locale-specific format;
//
////////////////////////////////////////////////////////////////////////////////
//
//  [char] chrToLower:  character to convert
//                  :  returns the converted character
//         
////////////////////////////////////////////////////////////////////////////////

char ConvertCharacterToLowerCase(char chrToLower)
{
	//  Locale, current
	locale locCurrent;
	//  The lower-case string to return
	char chrReturn = '\0';

	try
	{
		//  Convert character to lower-case
		chrReturn = tolower(chrToLower, locCurrent);
	}
	catch (exception ex)
	{
		cout << "ERROR [ConvertCharacterToLowerCase] Exception Code:  " << ex.what() << "\n";
	}

	return chrReturn;
}

////////////////////////////////////////////////////////////////////////////////
//
//  Converts a string to lower-case in a locale-specific format;
//
////////////////////////////////////////////////////////////////////////////////
//
//  [string] strToLower:  string to convert
//                         :  returns the converted string
//         
////////////////////////////////////////////////////////////////////////////////

string ConvertStringToLowerCase(string strToLower)
{
	//  The lower-case string to return
	string strReturn = "";

	try
	{
		//  Convert string to lower-case, character by character
		for (int nCount = 0; nCount < strToLower.length(); nCount++)
			strReturn += ConvertCharacterToLowerCase(strToLower[nCount]);
	}
	catch (exception ex)
	{
		cout << "ERROR [ConvertStringToLowerCase] Exception Code:  " << ex.what() << "\n";
	}

	return strReturn;
}

////////////////////////////////////////////////////////////////////////////////
//
//  Converts a long to a string
//         
////////////////////////////////////////////////////////////////////////////////
//
//  [long] lConvert:  long to convert
//                :  returns the string representation of the long
//         
////////////////////////////////////////////////////////////////////////////////

string ConvertLongToString(long lConvert)
{
	try
	{
		//  Convert
		#ifdef _WIN64
			return to_string(lConvert);
		#else
			//  String to return
			stringstream ssReturn;

			ssReturn << lConvert;

			return ssReturn.str();
		#endif
	}
	catch (exception ex)
	{
		cout << "ERROR [ConvertLongToString] Exception Code:  " << ex.what() << "\n";
	}

	return "";
}

////////////////////////////////////////////////////////////////////////////////
//
//  Converts a double to a string
//         
////////////////////////////////////////////////////////////////////////////////
//
//  [double] dConvert:  double to convert
//                  :  returns the string representation of the double
//         
////////////////////////////////////////////////////////////////////////////////

string ConvertDoubleToString(double dConvert)
{
	try
	{
		//  Convert
		#ifdef _WIN64
			return to_string(dConvert);
		#else
			//  String to return
			stringstream ssReturn;

			ssReturn << dConvert;

			return ssReturn.str();
		#endif
	}
	catch (exception ex)
	{
		cout << "ERROR [ConvertDoubleToString] Exception Code:  " << ex.what() << "\n";
	}

	return "";
}

////////////////////////////////////////////////////////////////////////////////
//
//  Converts a uint64_t to a string
//         
////////////////////////////////////////////////////////////////////////////////
//
//  [uint64_t] untConvert:  uint64_t to convert
//                      :  returns the string representation of the uint64_t
//         
////////////////////////////////////////////////////////////////////////////////

string ConvertUnsignedInt64ToString(uint64_t untConvert)
{
	try
	{
		//  Convert
		#ifdef _WIN64
			return to_string(untConvert);
		#else
			//  String to return
			stringstream ssReturn;

			ssReturn << untConvert;

			return ssReturn.str();
		#endif
	}
	catch (exception ex)
	{
		cout << "ERROR [ConvertUnsignedInt64ToString] Exception Code:  " << ex.what() << "\n";
	}

	return "";
}

////////////////////////////////////////////////////////////////////////////////
//
//  Converts a uint32_t to a string
//         
////////////////////////////////////////////////////////////////////////////////
//
//  [uint32_t] untConvert:  uint32_t to convert
//                      :  returns the string representation of the uint32_t
//         
////////////////////////////////////////////////////////////////////////////////

string ConvertUnsignedInt32ToString(uint32_t untConvert)
{
	try
	{
		//  Convert
		#ifdef _WIN64
			return to_string(untConvert);
		#else
			//  String to return
			stringstream ssReturn;

			ssReturn << untConvert;

			return ssReturn.str();
		#endif
	}
	catch (exception ex)
	{
		cout << "ERROR [ConvertUnsignedInt32ToString] Exception Code:  " << ex.what() << "\n";
	}

	return "";
}

////////////////////////////////////////////////////////////////////////////////
//
//  Converts a uint16_t to a string
//         
////////////////////////////////////////////////////////////////////////////////
//
//  [uint16_t] untConvert:  uint16_t to convert
//                      :  returns the string representation of the uint16_t
//        
////////////////////////////////////////////////////////////////////////////////

string ConvertUnsignedInt16ToString(uint16_t untConvert)
{
	try
	{
		//  Convert
		#ifdef _WIN64
			return to_string(untConvert);
		#else
			//  String to return
			stringstream ssReturn;

			ssReturn << untConvert;

			return ssReturn.str();
		#endif
	}
	catch (exception ex)
	{
		cout << "ERROR [ConvertUnsignedInt16ToString] Exception Code:  " << ex.what() << "\n";
	}

	return "";
}

////////////////////////////////////////////////////////////////////////////////
//
//  Converts an integer to a string
//         
////////////////////////////////////////////////////////////////////////////////
//
//  [int] nConvert:  integer to convert
//               :  returns the string representation of the integer
//         
////////////////////////////////////////////////////////////////////////////////

string ConvertIntToString(int nConvert)
{
	try
	{
		try
		{
			//  Convert
			#ifdef _WIN64
				return to_string(nConvert);
			#else
				//  String to return
				stringstream ssReturn;

				ssReturn << nConvert;

				return ssReturn.str();
			#endif
		}
		catch (exception ex)
		{
			cout << "ERROR [ConvertLongToString] Exception Code:  " << ex.what() << "\n";
		}

		return "";
	}
	catch (exception ex)
	{
		cout << "ERROR [ConvertIntToString] Exception Code:  " << ex.what() << "\n";
	}

	return "";
}

////////////////////////////////////////////////////////////////////////////////
//
//  Converts FNA to FA file format (first line is the header, second line contains all nucleotides) and extracts the FNA file accession; 
//      replaces '.' characters with '_' characters; replaces thymine (t) with uracil (u)
//         
////////////////////////////////////////////////////////////////////////////////
//
//  [string&] strFNAFileText:  FNA file text to convert
//  [string&] strFAFileText :  converted FA file text
//                         :  returns true, if successful; else, false
//         
////////////////////////////////////////////////////////////////////////////////

bool ConvertFNAtoFA(string& strFNAFileText, string& strFAFileText)
{
	//  Count the new-line-characters, max is one
	int nNewLineCount = 0;

	try
	{
		//  Clear the FA file text
		strFAFileText = "";

		//  If the file text is not empty, then convert the file
		if (!strFNAFileText.empty())
		{
			//  Strip all carriage-return characters, if any, and all new-line-characters except the first one delimiting the header line
			for (int nCount = 0; nCount < strFNAFileText.length(); nCount++)
			{
				//  Strip carriage-return characters, if any
				if (strFNAFileText[nCount] != '\r')
				{
					//  If new-line counter is zero OR the new-line counter is greater than zero and the current character is Not a new-line, concatenate the converted string
					if ((nNewLineCount == 0) || ((nNewLineCount > 0) && (strFNAFileText[nCount] != '\n')))
					{
						//  If not header, oncatenate the converted string, replacing thymine (t) with uracil (u)
						if ((nNewLineCount > 0) && (ConvertCharacterToLowerCase(strFNAFileText[nCount]) == 't'))
							strFAFileText += 'u';
						//  Concatenate the converted string
						else if (nNewLineCount > 0)
							strFAFileText += ConvertCharacterToLowerCase(strFNAFileText[nCount]);
						else
							strFAFileText += strFNAFileText[nCount];
					}

					//  If the current character is a new-line, then increment the new-line-character counter
					if ((nNewLineCount == 0) && (strFNAFileText[nCount] == '\n'))
						//  Increment the new-line-character counter
						nNewLineCount++;
				}
			}

			return true;
		}
		else
		{
			ReportTimeStamp("[ConvertFNAtoFA]", "ERROR:  Input FNA File Text is Empty");
		}
	}
	catch (exception ex)
	{
		cout << "ERROR [ConvertFNAtoFA] Exception Code:  " << ex.what() << "\n";
	}

	return false;
}

////////////////////////////////////////////////////////////////////////////////
//
//  Splits a string into a vector based on the provided delimiter (referenced for memory)
//         
////////////////////////////////////////////////////////////////////////////////
//
//  [const string &] strSplit :  string to split
//  [char] chrDelimiter       :  delimiter character
//  [vector<string> &] vTarget:  target vector<string>
//         
////////////////////////////////////////////////////////////////////////////////

void SplitString(const string &strSplit, char chrDelimiter, vector<string> &vTarget)
{
	try
	{
		//  The stringstream to process
		stringstream ssSplit(strSplit);
		//  The split string to deposit into the list<string>
		string strItem;

		//  Clear the list
		vTarget.clear();

		//  If the split string is not empty, then check delimiter value
		if (!strSplit.empty())
		{
			//  If the delimiter character is not empty, then split the string
			if (chrDelimiter != '\0')
			{
				//  Split the string into list<string>
				while (getline(ssSplit, strItem, chrDelimiter))
				{
					//  If the item is not empty
					if (!strItem.empty())
					{
						//  Add the string to the list<string>
						vTarget.push_back(strItem);
					}
				}
			}
			else
			{
				ReportTimeStamp("[SplitString]", "ERROR:  Delimiter character is Empty.");
			}
		}
		else
		{
			ReportTimeStamp("[SplitString]", "ERROR:  Input String is Empty");
		}
	}
catch (exception ex)
{
	cout << "ERROR [SplitString] Exception Code:  " << ex.what() << "\n";
}

return;
}

////////////////////////////////////////////////////////////////////////////////
//
//  Splits a string into a deque based on the provided delimiter (referenced for memory)
//         
////////////////////////////////////////////////////////////////////////////////
//
//  [const string &] strSplit :  string to split
//  [char] chrDelimiter       :  delimiter character
//  [deque<string> &] dTarget:  target vector<string>
//         
////////////////////////////////////////////////////////////////////////////////

void SplitString(const string &strSplit, char chrDelimiter, deque<string> &dTarget)
{
	try
	{
		//  The stringstream to process
		stringstream ssSplit(strSplit);
		//  The split string to deposit into the list<string>
		string strItem;

		//  Clear the list
		dTarget.clear();

		//  If the split string is not empty, then check delimiter value
		if (!strSplit.empty())
		{
			//  If the delimiter character is not empty, then split the string
			if (chrDelimiter != '\0')
			{
				//  Split the string into list<string>
				while (getline(ssSplit, strItem, chrDelimiter))
				{
					//  If the item is not empty
					if (!strItem.empty())
					{
						//  Add the string to the list<string>
						dTarget.push_back(strItem);
					}
				}
			}
			else
			{
				ReportTimeStamp("[SplitString]", "ERROR:  Delimiter character is Empty.");
			}
		}
		else
		{
			ReportTimeStamp("[SplitString]", "ERROR:  Input String is Empty");
		}
	}
	catch (exception ex)
	{
		cout << "ERROR [SplitString] Exception Code:  " << ex.what() << "\n";
	}

	return;
}

////////////////////////////////////////////////////////////////////////////////
//
//  Splits a string into a vector based on the provided delimiter (referenced for memory), leaves blank entries
//         
////////////////////////////////////////////////////////////////////////////////
//
//  [const string &] strSplit :  string to split
//  [char] chrDelimiter       :  delimiter character
//  [vector<string> &] vTarget:  target vector<string>
//         
////////////////////////////////////////////////////////////////////////////////

void SplitStringAllowEmptyEntries(const string &strSplit, string strDelimiter, vector<string> &vTarget)
{
	//  Character
	string strGet = "";
	//  Column entry
	string strEntry = "";

	try
	{
		//  Clear the list
		vTarget.clear();

		//  If the split string is not empty, then check delimiter value
		if (!strSplit.empty())
		{
			//  If the delimiter character is not empty, then split the string
			if (!strDelimiter.empty())
			{
				//  Iterate target string, push_back entries
				for (long lCount = 0; lCount <= strSplit.length(); lCount++)
				{
					//  Delimit entry
					if ((lCount == strSplit.length()) || (strSplit.substr(lCount, 1) == strDelimiter))
					{
						//  Add entry
						vTarget.push_back(strEntry);

						//  Clear entry
						strEntry = "";
					}
					else
						strEntry += strSplit.substr(lCount, 1);
				}
			}
			else
			{
				ReportTimeStamp("[SplitStringAllowEmptyEntries]", "ERROR:  Delimiter is Empty.");
			}
		}
		else
		{
			ReportTimeStamp("[SplitStringAllowEmptyEntries]", "ERROR:  Input String is Empty");
		}
	}
	catch (exception ex)
	{
		cout << "ERROR [SplitStringAllowEmptyEntries] Exception Code:  " << ex.what() << "\n";
	}

	return;
}

////////////////////////////////////////////////////////////////////////////////
//
//  Removes multiple file(s) in a list
//         
////////////////////////////////////////////////////////////////////////////////
//
//  [string] strInputFilePathName :  input file path name list
//  [string] strOutputFilePathName:  output list file path name
//  [string] strExtract              :  find text to replace
//  [string] strInsert           :  replace text upon find
//  [bool] bIsCaseSensitive       :  find is case sensitive, if true; else, not
//                               :  returns true, if successful; else, false
//         
////////////////////////////////////////////////////////////////////////////////

bool FileTextReplace(string strInputFilePathName, string strOutputFilePathName, string strExtract, string strInsert, bool bIsCaseSensitive)
{
	//  Input file text
	string strInputFileText = "";
	//  Output file text
	string strOutputFileText = "";

	try
	{
		//  If input file path name is not empty
		if (!strInputFilePathName.empty())
		{
			//  If output file path name is not empty
			if (!strOutputFilePathName.empty())
			{
				//  If output list file path name is not empty
				if (!strExtract.empty())
				{
					//  Ensure proper characters within find string
					strExtract = ReplaceInString(strExtract, "\\r", "\r", false);
					strExtract = ReplaceInString(strExtract, "\\n", "\n", false);
					strExtract = ReplaceInString(strExtract, "\\t", "\t", false);

					//  Ensure proper characters within replace string
					strInsert = ReplaceInString(strInsert, "\\r", "\r", false);
					strInsert = ReplaceInString(strInsert, "\\n", "\n", false);
					strInsert = ReplaceInString(strInsert, "\\t", "\t", false);

					//  Get input file text
					if (GetFileText(strInputFilePathName, strInputFileText))
					{
						//  If input file text is not empty
						if (!strInputFileText.empty())
						{
							//  Iterate string and replace
							for (long lCount = 0; lCount < strInputFileText.length(); lCount++)
							{
								//  If text matches find
								if (((!bIsCaseSensitive) && (strInputFileText.substr(lCount, strExtract.length()) == strExtract)) ||
									((bIsCaseSensitive) && (ConvertStringToLowerCase(strInputFileText.substr(lCount, strExtract.length())) == ConvertStringToLowerCase(strExtract))))
								{
									//  Replace text
									strOutputFileText += strInsert;

									//  Advance count
									lCount += (long)strExtract.length() - 1;
								}
								//  Else concatenate output
								else
									strOutputFileText += strInputFileText.substr(lCount, 1);
							}

							//  Write the output file
							if (!WriteFileText(strOutputFilePathName, strOutputFileText))
								ReportTimeStamp("[FileTextReplace]", "ERROR:  Output File [" + strOutputFilePathName + "] Write Failed");

							return true;
						}
						else
						{
							ReportTimeStamp("[FileTextReplace]", "ERROR:  Input File [" + strInputFilePathName + "] Text is Empty");
						}
					}
					else
					{
						ReportTimeStamp("[FileTextReplace]", "ERROR:  Input File [" + strInputFilePathName + "] Open Failed");
					}
				}
				else
				{
					ReportTimeStamp("[FileTextReplace]", "ERROR:  Find Text is Empty");
				}
			}
			else
			{
				ReportTimeStamp("[FileTextReplace]", "ERROR:  Output File Path Name is Empty");
			}
		}
		else
		{
			ReportTimeStamp("[FileTextReplace]", "ERROR:  Input File Path Name is Empty");
		}
	}
	catch (exception ex)
	{
		cout << "ERROR [FileTextReplace] Exception Code:  " << ex.what() << "\n";
	}

	return false;
}

////////////////////////////////////////////////////////////////////////////////
//
//  Differentiates two vectors via raw list comparison:
//
////////////////////////////////////////////////////////////////////////////////
//
//  [vector<string>&] vFrom :  from set
//  [vector<string>&] vTo   :  to set
//  [vector<string>&] vDiffs:  to set
//                          :  returns true, if successful; else, false
//         
////////////////////////////////////////////////////////////////////////////////

bool ListDiffVectors(vector<string>& vFrom, vector<string>& vTo, vector<string>& vDiffs)
{
	try
	{
		//  Clear output vector
		vDiffs.clear();

		//  Search from vector
		for (long lCountFrom = 0; lCountFrom < vFrom.size(); lCountFrom++)
		{
			//  Found if true; else false
			bool bFound = false;

			//  Search to vector
			for (long lCountTo = 0; lCountTo < vFrom.size(); lCountTo++)
			{
				//  If match, break
				if (vFrom[lCountFrom] == vTo[lCountTo])
				{
					break;
				}

				//  Increment count
				lCountTo++;
			}

			//  Add missing entry to output vector
			if (!bFound)
			{
				vDiffs.push_back(vFrom[lCountFrom]);
			}
		}

		return true;
	}
	catch (exception ex)
	{
		cout << "ERROR [ListDiffVectors] Exception Code:  " << ex.what() << "\n";
	}

	return false;
}

////////////////////////////////////////////////////////////////////////////////
//
//  Differentiates two vectors via set method:
//
////////////////////////////////////////////////////////////////////////////////
//
//  [vector<string>&] vFrom :  from set
//  [vector<string>&] vTo   :  to set
//  [vector<string>&] vDiffs:  to set
//                          :  returns true, if successful; else, false
//         
////////////////////////////////////////////////////////////////////////////////

bool VectorDiffVectors(vector<string>& vFrom, vector<string>& vTo, vector<string>& vDiffs)
{
	try
	{
		//  Sort vectors
		sort(vFrom.begin(), vFrom.end());
		sort(vTo.begin(), vTo.end());

		//  Diff vectors
		set_symmetric_difference(vFrom.begin(), vFrom.end(), vTo.begin(), vTo.end(), back_inserter(vDiffs));

		/*for (long lCount = vFrom.size() - 10; lCount < vFrom.size(); lCount++)
		{
			ReportTimeStamp("[VectorDiffVectors]", "NOTE:  From:  " + vFrom[lCount]);
		}

		for (long lCount = vTo.size() - 10; lCount < vTo.size(); lCount++)
		{
			ReportTimeStamp("[VectorDiffVectors]", "NOTE:  To:  " + vTo[lCount]);
		}

		for (long lCount = vDiffs.size() - 10; lCount < vDiffs.size(); lCount++)
		{
			ReportTimeStamp("[VectorDiffVectors]", "NOTE:  Diff:  " + vDiffs[lCount]);
		}*/

		return true;
	}
	catch (exception ex)
	{
		cout << "ERROR [VectorDiffVectors] Exception Code:  " << ex.what() << "\n";
	}

	return false;
}

////////////////////////////////////////////////////////////////////////////////
//
//  Differentiates two vectors via set method:
//
////////////////////////////////////////////////////////////////////////////////
//
//  [vector<string>&] vFrom :  from set
//  [vector<string>&] vTo   :  to set
//  [vector<string>&] vDiffs:  to set
//                          :  returns true, if successful; else, false
//         
////////////////////////////////////////////////////////////////////////////////

bool SetDiffVectors(vector<string>& vFrom, vector<string>& vTo, vector<string>& vDiffs)
{
	try
	{
		//  Sort vectors
		sort(vFrom.begin(), vFrom.end());
		sort(vTo.begin(), vTo.end());

		set<string> sFrom(vFrom.begin(), vFrom.end());
		set<string> sTo(vTo.begin(), vTo.end());

		//  Diff lists
		set_difference(sFrom.begin(), sFrom.end(), sTo.begin(), sTo.end(), back_inserter(vDiffs));

		/*for (long lCount = vFrom.size() - 10; lCount < vFrom.size(); lCount++)
		{
			ReportTimeStamp("[SetDiffVectors]", "NOTE:  From:  " + vFrom[lCount]);
		}

		for (long lCount = vTo.size() - 10; lCount < vTo.size(); lCount++)
		{
			ReportTimeStamp("[SetDiffVectors]", "NOTE:  To:  " + vTo[lCount]);
		}

		for (long lCount = vDiffs.size() - 10; lCount < vDiffs.size(); lCount++)
		{
			ReportTimeStamp("[SetDiffVectors]", "NOTE:  Diff:  " + vDiffs[lCount]);
		}*/

		return true;
	}
	catch (exception ex)
	{
		cout << "ERROR [SetDiffVectors] Exception Code:  " << ex.what() << "\n";
	}

	return false;
}

////////////////////////////////////////////////////////////////////////////////
//
//  Differentiates two vectors via multiset method:
//
////////////////////////////////////////////////////////////////////////////////
//
//  [vector<string>&] vFrom :  from set
//  [vector<string>&] vTo   :  to set
//  [vector<string>&] vDiffs:  to set
//                          :  returns diffed vector<string>, if successful; else, NULL
//         
////////////////////////////////////////////////////////////////////////////////

bool MultisetDiffVectors(vector<string>& vFrom, vector<string>& vTo, vector<string>& vDiffs)
{
	try
	{
		// Initialize hash tables
		unordered_multiset<string> umsFrom(vFrom.begin(), vFrom.end());
		unordered_multiset<string> umsTo(vTo.begin(), vTo.end());

		// Go through the first hash table
		for (auto autoFrom = umsFrom.cbegin(); autoFrom != umsFrom.cend();)
		{
			// Find the current item in the second hash table
			auto autoTo = umsTo.find(*autoFrom);

			// Is it present?
			if (autoTo != umsTo.end())
			{
				// If so, remove it from both hash tables
				autoFrom = umsFrom.erase(autoFrom);
				umsTo.erase(autoTo);
			}
			else
				++autoFrom;
		}

		// Create a vector of the union of the remaining items
		vDiffs = vector<string>(umsFrom.begin(), umsFrom.end());

		vDiffs.insert(vDiffs.end(), umsTo.begin(), umsTo.end());

		/*for (long lCount = vFrom.size() - 10; lCount < vFrom.size(); lCount++)
		{
			ReportTimeStamp("[MultisetDiffVectors]", "NOTE:  From:  " + vFrom[lCount]);
		}

		for (long lCount = vTo.size() - 10; lCount < vTo.size(); lCount++)
		{
			ReportTimeStamp("[MultisetDiffVectors]", "NOTE:  To:  " + vTo[lCount]);
		}

		for (long lCount = vDiffs.size() - 10; lCount < vDiffs.size(); lCount++)
		{
			ReportTimeStamp("[MultisetDiffVectors]", "NOTE:  Diff:  " + vDiffs[lCount]);
		}*/

		return true;
	}
	catch (exception ex)
	{
		cout << "ERROR [MultisetDiffVectors] Exception Code:  " << ex.what() << "\n";
	}

	return false;
}

////////////////////////////////////////////////////////////////////////////////
//
//  Gets the reverse compliment of a nucleotide sequence; assumes thymine (t) replaced with uracil (u)
//         
////////////////////////////////////////////////////////////////////////////////
//
//  [string] strIn       :  string to search
//  [string] strExtract  :  string to extract (replace)
//  [string] strInsert   :  string to insert
//  [bool] bCaseSensetive:  if true, regard case
//                      :  returns the replaced string
//         
////////////////////////////////////////////////////////////////////////////////

string ReplaceInString(string strIn, string strExtract, string strInsert, bool bCaseSensetive)
{
	//  The replaced string to return
	string strOut = "";

	try
	{
		//  If in string is not empty
		if (!strIn.empty())
		{
			//  If the extract string is not empty
			if (!strExtract.empty())
			{
				//  If in string contains at leats one replacement
				if (strIn.find_first_of(strIn) != string::npos)
				{
					//  Search and replace
					for (long lCount = 0; lCount < strIn.length(); lCount++)
					{
						//  If find text is found
						if (((bCaseSensetive) && (strIn.substr(lCount, strExtract.length()) == strExtract)) || (ConvertStringToLowerCase(strIn.substr(lCount, strExtract.length())) == ConvertStringToLowerCase(strExtract)))
						{
							//  Replace
							strOut += strInsert;

							//  Increment count
							lCount += (long)strExtract.length() - 1;
						}
						//  Concatenate out string
						else
							strOut += strIn[lCount];
					}
				}
			}
			else
			{
				ReportTimeStamp("[ReplaceInString]", "ERROR:  Search String is Empty.");
			}
		}
		else
		{
			ReportTimeStamp("[ReplaceInString]", "ERROR:  Input String is Empty");
		}
	}
	catch (exception ex)
	{
		cout << "ERROR [ReplaceInString] Exception Code:  " << ex.what() << "\n";
	}

	return strOut;
}

////////////////////////////////////////////////////////////////////////////////
//
//  Pads a string (prepend or append) with the indicated string
//         
////////////////////////////////////////////////////////////////////////////////
//
//  [string] strIn :  string to pad
//  [string] strPad:  pad string
//  [int] nLength  :  final output length
//  [bool] bPrepend:  prepend pad, if true; else, append
//         
////////////////////////////////////////////////////////////////////////////////

string PadString(string strIn, string strPad, int nLength, bool bPrepend)
{
	//  Padded string to return
	string strOut = "";

	try
	{
		if (!strPad.empty())
		{
			strOut = strIn;
			while(strOut.length() < nLength)
			{
				if (bPrepend)
					strOut = strPad + strOut;
				else
					strOut += strPad;
			}
		}
		else
		{
			ReportTimeStamp("[PadString]", "ERROR:  Pad String is Empty");
		}
	}
	catch (exception ex)
	{
		cout << "ERROR [PadString] Exception Code:  " << ex.what() << "\n";
	}

	return strOut;
}


////////////////////////////////////////////////////////////////////////////////
//
//  Scrubs an entry string of unwanted (known delimiter) character(s)
//         
////////////////////////////////////////////////////////////////////////////////
//
//  [string&] strIn :  string to pad
//         
////////////////////////////////////////////////////////////////////////////////

string ScrubDelimitedEntryString(string& strIn)
{
	//  Scrubbed string to return
	string strOut = "";
	//  Get string
	string strGet = "";

	try
	{
		//if (!strIn.empty())
		//{
			for (long lCount = 0; lCount < strIn.length(); lCount++)
			{
				strGet = strIn.substr(lCount, 1);

				if ((strGet != "~") && (strGet != "`") && (strGet != "^") && (strGet != "|"))
					strOut += strGet;
			}
		//}
		//else
		//{
		//	ReportTimeStamp("[ScrubDelimitedEntryString]", "ERROR:  Input String is Empty");
		//}
	}
	catch (exception ex)
	{
		cout << "ERROR [ScrubDelimitedEntryString] Exception Code:  " << ex.what() << "\n";
	}

	return strOut;
}

////////////////////////////////////////////////////////////////////////////////
//
//  Gets the reverse compliment of a nucleotide sequence; assumes thymine (t) replaced with uracil (u)
//         
////////////////////////////////////////////////////////////////////////////////
//
//  [string] strForward:  forward sequence to process
//         
////////////////////////////////////////////////////////////////////////////////

string ConvertToReverseCompliment(string strForward)
{
	//  The reverse compliment string to return
	string strReverse = "";
	//  The working character
	char chrGet = '\0';

	try
	{
		//  If the forward sequence is not empty, then process the sequence
		if (!strForward.empty())
		{
			//  Read the string backwards and concatenate the reverse compliment
			for (int nCount = (int)strForward.length() - 1; nCount >= 0; nCount--)
			{
				//  Get the current character
				chrGet = ConvertCharacterToLowerCase(strForward[nCount]);

				//  Concatenate the reverse compliment
				if (chrGet == 'a')
					strReverse += 'u';
				else if (chrGet == 't')
					strReverse += 'a';
				else if (chrGet == 'u')
					strReverse += 'a';
				else if (chrGet == 'g')
					strReverse += 'c';
				else if (chrGet == 'c')
					strReverse += 'g';
				else
					strReverse += 'n';
			}
		}
		else
		{
			ReportTimeStamp("[ConvertToReverseCompliment]", "ERROR:  Forward Sequence is Empty");
		}
	}
	catch (exception ex)
	{
		cout << "ERROR [ConvertToReverseCompliment] Exception Code:  " << ex.what() << "\n";
	}

	return strReverse;
}

////////////////////////////////////////////////////////////////////////////////
//
//  Gets the BIG accession from standardized BIG format file path name
//         
////////////////////////////////////////////////////////////////////////////////
//
//  [string] strFilePathName:  file path name to parse
//                          :  returns the accession
//         
////////////////////////////////////////////////////////////////////////////////

string GetAccessionFromBIGFilePathName(string strFilePathName)
{
	try
	{
		//  If the file path name is not empty
		if (!strFilePathName.empty())
		{
			//  Get base file name (no extensions)
			return GetBaseFileName(strFilePathName);
		}
		else
		{
			ReportTimeStamp("[GetAccessionFromBIGFilePathName]", "ERROR:  File Path Name is Empty");
		}
	}
	catch (exception ex)
	{
		cout << "ERROR [GetAccessionFromBIGFilePathName] Exception Code:  " << ex.what() << "\n";
	}

	return "";
}

////////////////////////////////////////////////////////////////////////////////
//
//  Gets the BIG accession from standardized BIG format file path name, extended for second file name extension
//         
////////////////////////////////////////////////////////////////////////////////
//
//  [string] strFilePathName:  file path name to parse
//  [int] nMaxExtensions    :  max extensions to include in extended accession
//                          :  returns the accession
//         
////////////////////////////////////////////////////////////////////////////////

string GetExtendedAccessionFromBIGFilePathName(string strFilePathName, int nMaxExtensions)
{
	//  File name to return
	string strExtendedAccession = "";
	//  Path delimiter character
	char chrPathDelimiter = strPathDelimiter;
	//  File path name parts
	vector<string> vFilePathNameParts;
	//  File name parts
	vector<string> vFileNameParts;

	try
	{
		//  If input file path name is not empty
		if (!strFilePathName.empty())
		{
			//  Split file path name parts
			SplitString(strFilePathName, chrPathDelimiter, vFilePathNameParts);

			//  If file path name parts is not empty
			if (vFilePathNameParts.size() >= 1)
			{
				//  Split file name parts
				SplitString(vFilePathNameParts[vFilePathNameParts.size() - 1], '.', vFileNameParts);

				//  If file name parts is not empty
				if (vFileNameParts.size() >= 1)
				{
					//  Iterate file path name parts
					for (int nCountFileParts = 0; ((nCountFileParts < (int)vFileNameParts.size() - 1) && (nCountFileParts <= nMaxExtensions)); nCountFileParts++)
					{
						//  Concatenate path
						if (!strExtendedAccession.empty())
							strExtendedAccession += "_";
						strExtendedAccession += vFileNameParts[nCountFileParts];
					}
				}
			}
		}
	}
	catch (exception ex)
	{
		cout << "ERROR [GetExtendedAccessionFromBIGFilePathName] Exception Code:  " << ex.what() << "\n";
	}

	return strExtendedAccession;
}

////////////////////////////////////////////////////////////////////////////////
//
//  Gets the FNA file accession from the FNA file header, replaces '.' characters with '_' characters
//         
////////////////////////////////////////////////////////////////////////////////
//
//  [string&] strFileText :  FNA file text to read
//                       :  returns the accession
//         
////////////////////////////////////////////////////////////////////////////////

string GetAccessionFromFileHeader(string& strFileText)
{
	//  The accession to return
	string strAccession = "";
	//  The header column count
	int nHeaderColumnCount = 0;
	//  The file
	vector<string> vFile;
	//  The header
	vector<string> vHeader;

	try
	{
		//  If the file text is not empty, then convert the file
		if (!strFileText.empty())
		{
			//  Split the file
			SplitString(strFileText, '\n', vFile);

			//  If file is set1
			if (vFile.size() > 0)
			{
				SplitString(vFile[0], '|', vHeader);

				//  If header is set
				if (vHeader.size() >= 4)
				{
					//  Set the accession
					strAccession = vHeader[3];

					//  Replace '.' with '_'
					strAccession = ReplaceInString(strAccession, ".", "_", false);
				}
				else
				{
					ReportTimeStamp("[GetAccessionFromFileHeader]", "ERROR:  Input FNA Header is Not set");
				}
			}
			else
			{
				ReportTimeStamp("[GetAccessionFromFileHeader]", "ERROR:  Input FNA File is Not Set");
			}
		}
		else
		{
			ReportTimeStamp("[GetAccessionFromFileHeader]", "ERROR:  Input FNA File Text is Empty");
		}
	}
	catch (exception ex)
	{
		cout << "ERROR [GetAccessionFromFileHeader] Exception Code:  " << ex.what() << "\n";
	}

	return strAccession;
}

////////////////////////////////////////////////////////////////////////////////
//
//  Gets the sequence from a BIG format .fa file
//         
////////////////////////////////////////////////////////////////////////////////
//
//  [string&] strFileText :  FA sequence file text to read
//                       :  returns the sequence
//         
////////////////////////////////////////////////////////////////////////////////

string GetSequenceFromFAFile(string& strFileText)
{
	//  The sequence to return
	string strSequence = "";
	//  The file
	vector<string> vFile;

	try
	{
		//  If the file text is not empty, then convert the file
		if (!strFileText.empty())
		{
			//  Split the file
			SplitString(strFileText, '\n', vFile);

			//  If file is set1
			if (vFile.size() > 1)
			{
				strSequence = vFile[1];
			}
			else
			{
				ReportTimeStamp("[GetSequenceFromFAFile]", "ERROR:  Input FA File is Not Properly Formatted");
			}
		}
		else
		{
			ReportTimeStamp("[GetSequenceFromFAFile]", "ERROR:  Input FA File Text is Empty");
		}
	}
	catch (exception ex)
	{
		cout << "ERROR [GetSequenceFromFAFile] Exception Code:  " << ex.what() << "\n";
	}

	return strSequence;
}

////////////////////////////////////////////////////////////////////////////////
//
//  Gets the ORF length a BIG format .ORF file path name
//         
////////////////////////////////////////////////////////////////////////////////
//
//  [string] strORFFilePathName:  ORF file path name to parse
//                             :  returns the ORF length, if successful; else, -1
//         
////////////////////////////////////////////////////////////////////////////////

int GetLengthFromORFFilePathName(string strORFFilePathName)
{
	//  File name parts
	vector<string> vFileNameParts;
	//  ORF parts
	vector<string> vORFParts;

	try
	{
		//  If file path name is not empty
		if (!strORFFilePathName.empty())
		{
			// Split string, get file path name parts
			SplitString(GetFileName(strORFFilePathName), '.', vFileNameParts);

			//  If file name parts is not empty
			if (vFileNameParts.size() >= 2)
			{
				// Split string, get file name parts
				SplitString(vFileNameParts[1], '_', vORFParts);

				//  If ORF parts is not empty
				if (vORFParts.size() == 3)
				{
					// strORFFilePathName = strBasePathName + "." + strDirection + strFrameNumber + "_" + strStart + "_" + strStop + ".ORF";

					//  ORF start
					int nStart = 0;
					//  ORF stop
					int nStop = 0;

					//  Get ORF start
					stringstream(vORFParts[1]) >> nStart;

					//  Get ORF stop
					stringstream(vORFParts[2]) >> nStop;

					return nStop - nStart;
				}
			}
		}
		else
		{
			ReportTimeStamp("[GetLengthFromORFFilePathName]", "ERROR:  .ORF File Path Name is Empty");
		}
	}
	catch (exception ex)
	{
		cout << "ERROR [GetLengthFromORFFilePathName] Exception Code:  " << ex.what() << "\n";
	}

	return -1;
}

////////////////////////////////////////////////////////////////////////////////
//
//  Parses a BIG .ORF file path name into ORF sequence, orientation, start, stop and length values
//         
////////////////////////////////////////////////////////////////////////////////
//
//  [string] strORFFilePathName:  ORF file path name to parse
//  [string&] strOrientation   :  ORF frame orientation and number to return
//  [long&] lStart             :  ORF start position to return
//  [long&] lStop              :  ORF stop position to return
//  [long&] lLength            :  ORF length to return
//                             :  returns true, if successful; else, false
//         
////////////////////////////////////////////////////////////////////////////////

bool ParseORFFilePathName(string strORFFilePathName, string& strOrientation, long& lStart, long& lStop, long& lLength)
{
	//  File name parts
	vector<string> vFileNameParts;
	//  ORF parts
	vector<string> vORFParts;
	//  FA file text
	string strFAFileText = "";

	try
	{
		//  If file path name is not empty
		if (!strORFFilePathName.empty())
		{
			// Split string, get file path name parts
			SplitString(GetFileName(strORFFilePathName), '.', vFileNameParts);

			//  If file name parts is not empty
			if (vFileNameParts.size() >= 2)
			{
				// Split string, get file name parts
				SplitString(vFileNameParts[1], '_', vORFParts);

				//  If ORF parts is not empty
				if (vORFParts.size() == 3)
				{
					// strORFFilePathName = strBasePathName + "." + strDirection + strFrameNumber + "_" + strStart + "_" + strStop + ".ORF";

					//  Get orientation
					strOrientation = vORFParts[0];

					//  Get ORF start
					stringstream(vORFParts[1]) >> lStart;

					//  Get ORF stop
					stringstream(vORFParts[2]) >> lStop;

					//  Calculate length
					lLength = lStop - lStart;

					return true;
				}
			}
		}
	}
	catch (exception ex)
	{
		cout << "ERROR [ParseORFFilePathName] Exception Code:  " << ex.what() << "\n";
	}

	return false;
}

////////////////////////////////////////////////////////////////////////////////
//
//  Returns the message associated with a given error number
//         
////////////////////////////////////////////////////////////////////////////////
//
//  [string] strORFFilePathName:  .ORF file path name to convert
//                             :  returns the .fa file path name, if successful; else, an empty string
//         
////////////////////////////////////////////////////////////////////////////////

string ConvertORFToFAFilePathName(string strORFFilePathName)
{
	//  .fa file path name to return
	string strFAFilePathName = "";

	try
	{
		//  File path name parts
		vector<string> vFilePathNameParts;
		//  Delmiter
		char chrDelimiter = strPathDelimiter;

		//  If .ORF file path name is not empty
		if (!strORFFilePathName.empty())
		{
			//  Split file path name
			SplitString(strORFFilePathName, chrDelimiter, vFilePathNameParts);

			//  If file path name parts exist
			if (vFilePathNameParts.size() > 0)
			{
				//  Concatenate .fa file path name
				for (int nCount = 0; nCount < (int)vFilePathNameParts.size(); nCount++)
				{
					if ((ConvertStringToLowerCase(vFilePathNameParts[nCount]) == "vrna") || (ConvertStringToLowerCase(vFilePathNameParts[nCount]) == "rosa") || (ConvertStringToLowerCase(vFilePathNameParts[nCount]) == "temp"))
					{
						if((nCount-1 >= 0) && (nCount-1 < vFilePathNameParts.size()))
							strFAFilePathName += chrDelimiter + vFilePathNameParts[nCount-1];

						nCount = (int)vFilePathNameParts.size();
					}
					else
					{
						strFAFilePathName += chrDelimiter + vFilePathNameParts[nCount];
					}
				}

				//  Append extension
				strFAFilePathName += ".fa";
			}
			else
			{
				ReportTimeStamp("[ConvertORFToFAFilePathName]", "ERROR:  .ORF File Path Name [" + strORFFilePathName + "] Parse Failed");
			}
		}
		else
		{
			ReportTimeStamp("[ConvertORFToFAFilePathName]", "ERROR:  .ORF File Path Name is Empty");
		}
	}
	catch (exception ex)
	{
		cout << "ERROR [ConvertORFToFAFilePathName] Exception Code:  " << ex.what() << "\n";
	}

	return strFAFilePathName;
}

////////////////////////////////////////////////////////////////////////////////
//
//  Returns the message associated with a given error number
//         
////////////////////////////////////////////////////////////////////////////////
//
//  [int] nError:  error number to interpret
//              :  returns the error message, if successful; else, returns "unknown error"
//         
////////////////////////////////////////////////////////////////////////////////

string GetErrorMessage(int nError)
{
	try
	{
		if (nError == E2BIG)
			return "The total number of bytes in the environment(envp) and argument list(argv) is too large.";
		else if (nError == EACCES)
			return "Search permission is denied on a component of the path prefix of filename or the name of a script interpreter. (See also path_resolution(7).);\n-OR-\nThe file or a script interpreter is not a regular file.\n-OR-\nExecute permission is denied for the file or a script or ELF interpreter.\n-OR\nThe file system is mounted noexec.";
		else if (nError == EFAULT)
			return "Filename points outside your accessible address space.";
		else if (nError == EINVAL)
			return "An ELF executable had more than one PT_INTERP segment(i.e., tried to name more than one interpreter).";
		else if (nError == EIO)
			return "n I / O error occurred.";
		else if (nError == EISDIR)
			return "An ELF interpreter was a directory.";
		//else if (nError == ELIBBAD)
		//	return "An ELF interpreter was not in a recognized format.";
		else if (nError == ELOOP)
			return "Too many symbolic links were encountered in resolving filename or the name of a script or ELF interpreter.";
		else if (nError == EMFILE)
			return "The process has the maximum number of files open.";
		else if (nError == ENAMETOOLONG)
			return "filename is too long.";
		else if (nError == ENFILE)
			return "The system limit on the total number of open files has been reached.";
		else if (nError == ENOENT)
			return "The file filename or a script or ELF interpreter does not exist, or a shared library needed for file or interpreter cannot be found.";
		else if (nError == ENOEXEC)
			return "An executable is not in a recognized format, is for the wrong architecture, or has some other format error that means it cannot be executed.";
		if (nError == ENOMEM)
			return "Insufficient kernel memory was available.";
		else if (nError == ENOTDIR)
			return "A component of the path prefix of filename or a script or ELF interpreter is not a directory.";
		else if (nError == EPERM)
			return "The file system is mounted nosuid, the user is not the superuser, and the file has the set - user - ID or set - group - ID bit set.";
		else if (nError == EPERM)
			return "The process is being traced, the user is not the superuser and the file has the set - user - ID or set - group - ID bit set.";
		else if (nError == ETXTBSY)
			return "Executable was open for writing by one or more processes.";
		else
			return "Unknown error.";
	}
	catch (exception ex)
	{
		cout << "ERROR [GetErrorMessage] Exception Code:  " << ex.what() << "\n";
	}

	return "";
}
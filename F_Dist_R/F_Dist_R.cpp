// F_Dist_R.cpp : Defines the entry point for the console application.
//

#include "F_Dist_R.h"
#include "Base123_Utilities.h"
#include "Base123_Catalog_Entry.h"
#include "Base123_Catalog.h"
#include "Base123_FDistance.h"

using namespace std;

#include <sys/stat.h>
#include <sstream>

////////////////////////////////////////////////////////////////////////////////
//
//  See ReportBase123Help() function for algorithm details;
//         
////////////////////////////////////////////////////////////////////////////////
//
//  [int] nArgumentCount:  argument count, space delimited
//  [char*] chpArguments:  arguments, space cropped
//                     :  returns the exit or return code
//         
////////////////////////////////////////////////////////////////////////////////

int main(int nArgumentCount, char* chpArguments[])
{
	//  Parameter report
	string strReport = "";
	//  Arguments vector
	vector<string> vArgs;

	try
	{
		#ifdef _WIN64
		#else
			//  Set the umask
			umask(0007);
		#endif

		// Command line must contain at least (list file and output folder [or use base path as output folder])
		//    or it must call for help, only
		if (nArgumentCount < 2)
		{
			// Update console; suggest help
			ReportTimeStamp(chpArguments[0], "ERROR: Use the Following Command Switch for Assistance: -help [-h]");

			return -1;
		}
		// Process switches and parameters
		else
		{
			//  Replace empty arguments with empty strings
			for (int nCount = 0; nCount < nArgumentCount; nCount++)
			{
				if (ConvertStringToLowerCase(chpArguments[nCount]) == "^mt^")
				{
					string strMT = "";
					vArgs.push_back(strMT);
				}
				else
					vArgs.push_back(chpArguments[nCount]);
			}

			//  Report  help
			if ((ConvertStringToLowerCase(vArgs[1]) == "-help") || (ConvertStringToLowerCase(vArgs[1]) == "-h"))
			{
				//ReportBase123Help();

				return 0;
			}


			//  Update console; start application;
			ReportTimeStamp(vArgs[0], "F_Dist_R Start: ");

			//  Report parameters
			strReport = "[" + ConvertIntToString(nArgumentCount) + "] ";
			for (int nCount = 0; nCount < nArgumentCount; nCount++)
			{
				strReport += vArgs[nCount];
				strReport += "; ";
			}
			strReport += "\n";

			//  Report to console
			ReportTimeStamp("Parameters:  ", strReport);

			//  Test
			if ((ConvertStringToLowerCase(vArgs[1]) == "-test") || (ConvertStringToLowerCase(vArgs[1]) == "-t"))
			{
			}
			//  Perform F-Distance analysis
			else if ((ConvertStringToLowerCase(vArgs[1]) == "-perform_fdistance_analysis") || (ConvertStringToLowerCase(vArgs[1]) == "-pfda"))
			{
				//  Usage is Base123 <switch> <arg1> <arg2> <arg3> <arg4> <arg5> <arg6> <arg7> <arg8> <arg9> <arg10> <arg11> <arg12> <arg13> <arg14> <arg15> <arg16> <arg17> <arg18>
				//    -perform_fdistance_analysis [-pfda]
				//         <output_table_file_path_name>
				//         <background_input_file_path_name_list>
				//         <background_input_file_path_name_transform>
				//         <background_catalog_file_path_name>
				//         <maximum_background_catalog_size>
				//         -background_unidirect [-bu]
				//              ...OR
				//                   -background_bidirect [-bb]
				//         -background_allow_unknowns [-bau]
				//              ...OR
				//                   -background_disallow_unknowns [-bdu]
				//         <background_error_file_path_name>
				//         <foreground_input_file_path_name_list>
				//         <foreground_input_file_path_name_transform>
				//         <foreground_catalog_file_path_name>
				//         <maximum_background_catalog_size>
				//         -foreground_unidirect [-fu]
				//              ...OR
				//                   -foreground_bidirect [-fb]
				//         <output_file_name_suffix>
				//         <foreground_error_file_path_name>
				//         <nmer_length>
				//         <max_processors>

				if (nArgumentCount == 20)
				{
					bool bResult = false;
					string strOutputTableFilePathName = "";
					string strBackgroundFilePathNameList = "";
					string strBackgroundInputFilePathNameTransform = "";
					string strBackgroundCatalogFilePathName = "";
					long lMaxBackgroundCatalogSize = 0;
					bool bBackgroundBidirect = false;
					bool bBackgroundAllowUnknowns = false;
					string strBackgroundErrorFilePathName = "";
					string strForegroundFilePathNameList = "";
					string strForegroundInputFilePathNameTransform = "";
					string strForegroundCatalogFilePathName = "";
					long lMaxForegroundCatalogSize = 0;
					bool bForegroundBidirect = false;
					bool bForegroundAllowUnknowns = false;
					string strForegroundErrorFilePathName = "";
					string strOutputFileNameSuffix = "";
					int nNMerLength = 0;
					int nMaxProcs = 0;

					strOutputTableFilePathName = vArgs[2];
					strBackgroundFilePathNameList = vArgs[3];
					strBackgroundInputFilePathNameTransform = vArgs[4];
					strBackgroundCatalogFilePathName = vArgs[5];
					stringstream(vArgs[6]) >> lMaxBackgroundCatalogSize;
					if ((ConvertStringToLowerCase(vArgs[7]) == "-background_bidirect") || (ConvertStringToLowerCase(vArgs[7]) == "-bb"))
						bBackgroundBidirect = true;
					if ((ConvertStringToLowerCase(vArgs[8]) == "-background_allow_unknowns") || (ConvertStringToLowerCase(vArgs[8]) == "-bau"))
						bBackgroundAllowUnknowns = true;
					strBackgroundErrorFilePathName = vArgs[9];
					strForegroundFilePathNameList = vArgs[10];
					strForegroundInputFilePathNameTransform = vArgs[11];
					strForegroundCatalogFilePathName = vArgs[12];
					stringstream(vArgs[13]) >> lMaxForegroundCatalogSize;
					if ((ConvertStringToLowerCase(vArgs[14]) == "-foreground_bidirect") || (ConvertStringToLowerCase(vArgs[14]) == "-fb"))
						bForegroundBidirect = true;
					if ((ConvertStringToLowerCase(vArgs[15]) == "-foreground_allow_unknowns") || (ConvertStringToLowerCase(vArgs[15]) == "-fau"))
						bForegroundAllowUnknowns = true;
					strOutputFileNameSuffix = vArgs[16];
					strForegroundErrorFilePathName = vArgs[17];
					stringstream(vArgs[18]) >> nNMerLength;
					stringstream(vArgs[19]) >> nMaxProcs;

					if ((nNMerLength == 8) || (nNMerLength == 16))
					{
						if (PerformFDistanceAnalysis(strOutputTableFilePathName, strBackgroundFilePathNameList, strBackgroundInputFilePathNameTransform,
							strBackgroundCatalogFilePathName, lMaxBackgroundCatalogSize, bBackgroundBidirect, bBackgroundAllowUnknowns, strBackgroundErrorFilePathName,
							strForegroundFilePathNameList, strForegroundInputFilePathNameTransform, strForegroundCatalogFilePathName, lMaxForegroundCatalogSize,
							bForegroundBidirect, bForegroundAllowUnknowns, strOutputFileNameSuffix, strForegroundErrorFilePathName, nNMerLength, nMaxProcs))
						{
							ReportTimeStamp(vArgs[0], "ERROR:  F-Distance Analysis Failed");

							return -1;
						}
					}
					else
					{
						ReportTimeStamp(vArgs[0], "ERROR:  F-Distance nMer Length Must be [8, 16]:  Use -help [-h] Switch for Assistance");

						return -1;
					}
				}
				else
				{
					ReportTimeStamp(vArgs[0], "ERROR:  Command Line is Not Properly Formatted to Perform F-Distance Analysis:  Use -help [-h] Switch for Assistance");

					return -1;
				}
			}
			else
			{
				//  Report unrecognized switch set
				string strMsg("ERROR:  Invalid Command Switch [^ARG^]");
				strMsg.replace(strMsg.find("^ARG^"), 5, vArgs[1]);
				ReportTimeStamp(vArgs[0], (char*)strMsg.c_str());

				return -1;
			}

			// Update console; end application;
			ReportTimeStamp(vArgs[0], "F_Dist_R Finish");
		}
	}
	catch (exception ex)
	{
		cout << "ERROR [F_Dist_R] Exception Code (Check Parameter Arguments/Counts - Index Overflow?):  " << ex.what() << "\n";

		return -1;
	}

	return 0;
}

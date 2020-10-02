// Base123_FDistance.cpp : Performs Base123 Fofanov Distance Genomic Analysis

////////////////////////////////////////////////////////////////////////////////
//
//  Performs Base123 Fofanov Distance Genomic Analysis; compares a foreground genome to a background genome
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

#include "F_Dist_R.h"
#include "Base123_Catalog_Entry.h"
#include "Base123_Catalog.h"
#include "Base123_Utilities.h"
#include "Base123_FDistance_16.h"
#include "Base123_FDistance_32.h"
#include "Base123_FDistance.h"

#include <math.h>
#include <sstream>
#include <omp.h>

////////////////////////////////////////////////////////////////////////////////
//
//  Filters a BIG FA Format File for F-Distance Suitability
//         
////////////////////////////////////////////////////////////////////////////////
//
//  [string] strInputFilePathName  :  input FA format file path name
//  [bool] bUseStrictFilter        :  use strict filter, if true; else, use lax filter
//  [int] nMaxPolyLimit            :  maximum poly-character limit
//  [string&] strAcceptListFileText:  output list file text
//  [string&] strRejectListFileText:  output list file text
//                                :  returns true, if successful; else, false
//         
////////////////////////////////////////////////////////////////////////////////

bool FilterFileForFDistance(string strInputFilePathName, bool bUseStrictFilter, int nMaxPolyLimit, string& strAcceptListFileText, string& strRejectListFileText)
{
	//  Input file text
	string strInputFileText = "";
	//  File header and text
	vector<string> vFileParts;
	//  Reject entry
	string strRejectEntry = "";
	//  Poly string
	string strPoly = "";
	//  Accepted, if true; else, rejected
	bool bAccept = true;

	try
	{
		//  If input file path name is not empty
		if (!strInputFilePathName.empty())
		{
			//  Get file text
			if (GetFileText(strInputFilePathName, strInputFileText))
			{
				//  Get file parts
				SplitString(strInputFileText, '\n', vFileParts);

				//  If  the file parts are not empty
				if ((vFileParts.size() >= 2) && (!vFileParts[1].empty()))
				{
					//  Set reject entry leader
					string strRejectEntry = strInputFilePathName;

					/*
					R	A or G	puRine
					Y	C, T or U	pYrimidines
					K	G, T or U	bases which are Ketones
					M	A or C	bases with aMino groups
					S	C or G	Strong interaction
					W	A, T or U	Weak interaction
					B	not A (i.e. C, G, T or U)	B comes after A
					D	not C (i.e. A, G, T or U)	D comes after C
					H	not G (i.e., A, C, T or U)	H comes after G
					V	neither T nor U (i.e. A, C or G)	V comes after U
					N	A C G T U	Nucleic acid
					-	gap of indeterminate length
					*/

					// "File~Gap (y/n)~Poly-n (y/n)~
					// "Strict - n(y / n)~Strict - r(y / n)~Strict - y(y / n)~Strict - k(y / n)~Strict - m(y / n)~Strict - s(y / n)~Strict - w(y / n)~Strict - b(y / n)~Strict - d(y / n)~Strict - h(y / n)~Strict - v(y / n)\n";

					if (vFileParts[1].find_first_of("-") != string::npos)
					{
						bAccept = false;
						strRejectEntry += "~y";
					}
					else
						strRejectEntry += "~n";

					strPoly = PadString("", "n", nMaxPolyLimit, false);
					if (vFileParts[1].find(strPoly) != string::npos)
					{
						bAccept = false;
						strRejectEntry += "~y";
					}
					else
						strRejectEntry += "~n";

					if (bUseStrictFilter)
					{
						if (vFileParts[1].find_first_of("n") != string::npos)
						{
							bAccept = false;
							strRejectEntry += "~y";
						}
						else
							strRejectEntry += "~n";

						if (vFileParts[1].find_first_of("r") != string::npos)
						{
							bAccept = false;
							strRejectEntry += "~y";
						}
						else
							strRejectEntry += "~n";

						if (vFileParts[1].find_first_of("y") != string::npos)
						{
							bAccept = false;
							strRejectEntry += "~y";
						}
						else
							strRejectEntry += "~n";

						if (vFileParts[1].find_first_of("k") != string::npos)
						{
							bAccept = false;
							strRejectEntry += "~y";
						}
						else
							strRejectEntry += "~n";

						if (vFileParts[1].find_first_of("m") != string::npos)
						{
							bAccept = false;
							strRejectEntry += "~y";
						}
						else
							strRejectEntry += "~n";

						if (vFileParts[1].find_first_of("s") != string::npos)
						{
							bAccept = false;
							strRejectEntry += "~y";
						}
						else
							strRejectEntry += "~n";

						if (vFileParts[1].find_first_of("w") != string::npos)
						{
							bAccept = false;
							strRejectEntry += "~y";
						}
						else
							strRejectEntry += "~n";

						if (vFileParts[1].find_first_of("b") != string::npos)
						{
							bAccept = false;
							strRejectEntry += "~y";
						}
						else
							strRejectEntry += "~n";

						if (vFileParts[1].find_first_of("d") != string::npos)
						{
							bAccept = false;
							strRejectEntry += "~y";
						}
						else
							strRejectEntry += "~n";

						if (vFileParts[1].find_first_of("h") != string::npos)
						{
							bAccept = false;
							strRejectEntry += "~y";
						}
						else
							strRejectEntry += "~n";

						if (vFileParts[1].find_first_of("v") != string::npos)
						{
							bAccept = false;
							strRejectEntry += "~y";
						}
						else
							strRejectEntry += "~n";
					}

					if (bAccept)
						strAcceptListFileText = strInputFilePathName;
					else
						strRejectListFileText = strRejectEntry;

					return true;
				}
				else
				{
					ReportTimeStamp("[FilterFileForFDistance]", "ERROR:  Input File [" + strInputFilePathName + "] is Not Properly Formatted");
				}
			}
			else
			{
				ReportTimeStamp("[FilterFileForFDistance]", "ERROR:  Input File [" + strInputFilePathName + "] Open Failed");
			}
		}
		else
		{
			ReportTimeStamp("[FilterFileForFDistance]", "ERROR:  Input List File Path Name is Empty");
		}
	}
	catch (exception ex)
	{
		cout << "ERROR [FilterFileForFDistance] Exception Code:  " << ex.what() << "\n";
	}

	return false;
}

////////////////////////////////////////////////////////////////////////////////
//
//  Filters a BIG FA Format File List for F-Distance Suitability
//         
////////////////////////////////////////////////////////////////////////////////
//
//  [string] strInputFilePathNameList       :  input file path name list to process
//  [string] strInputFilePathNameTransform  :  input file path name transform (includes string replacements, see help)
//  [bool] bUseStrictFilter                 :  use strict filter, if true; else, use lax filter
//  [int] nMaxPolyLimit                     :  maximum poly-character limit
//  [string] strAcceptOutputListFilePathName:  accept output list file path name
//  [string] strRejectOutputListFilePathName:  reject output list file path name
//  [int] nMaxProcs                         :  maximum processors for openMP
//                                         :  returns true, if successful; else, false
//         
////////////////////////////////////////////////////////////////////////////////

bool ListFilterForFDistance(string strInputFilePathNameList, string strInputFilePathNameTransform, bool bUseStrictFilter, int nMaxPolyLimit, string strAcceptOutputListFilePathName,
	string strRejectOutputListFilePathName, int nMaxProcs)
{
	//  List lock
	omp_lock_t lockList;
	//  List file text
	string strInputListFileText = "";
	//  File path names extracted from the file path name list file
	vector<string> vFilePathNames;
	//  Reject file header
	string strHeader = "";
	//  Acceptions
	vector<string> vAcceptListEntries;
	//  Accept file text
	string strAcceptListFileText = "";
	//  Rejections
	vector<string> vRejectListEntries;
	//  Reject file text
	string strRejectListFileText = "";

	try
	{
		//  If input file path name list is not empty
		if (!strInputFilePathNameList.empty())
		{
			//  Get list file text
			if (GetFileText(strInputFilePathNameList, strInputListFileText))
			{
				//  Split file text into file list
				SplitString(strInputListFileText, '\n', vFilePathNames);

				//  If file path names exist
				if (vFilePathNames.size() > 0)
				{
					//  Initialize accept/reject containers
					vAcceptListEntries.resize(vFilePathNames.size());
					vRejectListEntries.resize(vFilePathNames.size());

					//  Initialize time stamp lock
					omp_init_lock(&lockList);

					//  Declare team size
					#pragma omp parallel num_threads(nMaxProcs)
					{
						#pragma omp for
						for (long lCount = 0; lCount < vFilePathNames.size(); lCount++)
						{
							//  Test max procs
							if (lCount == 0)
							{
								omp_set_lock(&lockList);
								ReportTimeStamp("[ListFilterForFDistance]", "NOTE:  Thread Count = " + ConvertIntToString(omp_get_num_threads()));
								omp_unset_lock(&lockList);
							}

							//  Update for timestamp every 10,000 files
							if (lCount % 10000 == 0)
							{
								omp_set_lock(&lockList);
								ReportTimeStamp("[ListFilterForFDistance]", "NOTE:  Processing Entry [" + ConvertLongToString(lCount) + "] [" + vFilePathNames[lCount] + "]");
								omp_unset_lock(&lockList);
							}

							//  If the file path name is not empty, then proceed
							if (!vFilePathNames[lCount].empty())
							{
								//  Working file path name
								string strWorkingFilePathName = "";

								//  If input file path name transform is not empty
								if (!strInputFilePathNameTransform.empty())
									strWorkingFilePathName = TransformFilePathName(vFilePathNames[lCount], strInputFilePathNameTransform, "");
								else
									strWorkingFilePathName = vFilePathNames[lCount];

								//  Filter file for F-Distance suitability
								if (!FilterFileForFDistance(strWorkingFilePathName, bUseStrictFilter, nMaxPolyLimit, vAcceptListEntries[lCount], vRejectListEntries[lCount]))
								{
									vRejectListEntries[lCount] += strWorkingFilePathName;

									omp_set_lock(&lockList);
									ReportTimeStamp("[ListFilterForFDistance]", "ERROR:  Input File [" + strWorkingFilePathName + "] Filter Failed");
									omp_unset_lock(&lockList);
								}
							}
						}
					}

					//  Destroy time stamp lock
					omp_destroy_lock(&lockList);

					//  Concatenate accepted/rejected file list file text
					for (long lCount = 0; lCount < vFilePathNames.size(); lCount++)
					{
						if (!vAcceptListEntries[lCount].empty())
							strAcceptListFileText += GetFileName(vAcceptListEntries[lCount]) + "\n";

						if (!vRejectListEntries[lCount].empty())
							strRejectListFileText += GetFileName(vRejectListEntries[lCount]) + "\n";
					}

					//  If accept file text is not empty
					if (!strAcceptOutputListFilePathName.empty())
					{
						//  Write the file
						if (!WriteFileText(strAcceptOutputListFilePathName, strAcceptListFileText))
							ReportTimeStamp("[ListFilterForFDistance]", "ERROR:  Accepted List File [" + strAcceptOutputListFilePathName + "] Write Failed");
					}

					//  Write rejected file list file text
					if (!strRejectOutputListFilePathName.empty())
					{
						/*
						R	A or G	puRine
						Y	C, T or U	pYrimidines
						K	G, T or U	bases which are Ketones
						M	A or C	bases with aMino groups
						S	C or G	Strong interaction
						W	A, T or U	Weak interaction
						B	not A (i.e. C, G, T or U)	B comes after A
						D	not C (i.e. A, G, T or U)	D comes after C
						H	not G (i.e., A, C, T or U)	H comes after G
						V	neither T nor U (i.e. A, C or G)	V comes after U
						N	A C G T U	Nucleic acid
						-	gap of indeterminate length
						*/

						strHeader = "File~Gap (y/n)~Poly-n (y/n)~";
						strHeader += "Strict - n(y / n)~Strict - r(y / n)~Strict - y(y / n)~Strict - k(y / n)~Strict - m(y / n)~Strict - s(y / n)~Strict - w(y / n)~Strict - b(y / n)~Strict - d(y / n)~Strict - h(y / n)~Strict - v(y / n)\n";
						strRejectListFileText = strHeader + strRejectListFileText;

						//  Write the file
						if (!WriteFileText(strRejectOutputListFilePathName, strRejectListFileText))
							ReportTimeStamp("[ListFilterForFDistance]", "ERROR:  Rejected List File [" + strRejectOutputListFilePathName + "] Write Failed");
					}

					return true;
				}
				else
				{
					ReportTimeStamp("[ListFilterForFDistance]", "ERROR:  Input List File [" + strInputFilePathNameList + "] File Text is Empty");
				}
			}
			else
			{
				ReportTimeStamp("[ListFilterForFDistance]", "ERROR:  Input List File [" + strInputFilePathNameList + "] Open Failed");
			}
		}
		else
		{
			ReportTimeStamp("[ListFilterForFDistance]", "ERROR:  Input List File Path Name is Empty");
		}
	}
	catch (exception ex)
	{
		cout << "ERROR [ListFilterForFDistance] Exception Code:  " << ex.what() << "\n";
	}

	return false;
}

////////////////////////////////////////////////////////////////////////////////
//
//  Tabulates a single F-Distance output file and shuffle-file set
//         
////////////////////////////////////////////////////////////////////////////////
//
//  [string&] strInputFilePathName   :  BIG format .fa file path name
//  [string&] strAccession           :  NCBI accession number
//  [int] nOutputCount               :  original Shuffler output file count
//  [string&] strOutputFileNameSuffix:  original F-Distance output file name suffix
//  [string&] strTableEntry          :  table entry text
//  [string&] strErrorEntry          :  error entry text
//                                   :  returns true, if successful; else, false
//         
////////////////////////////////////////////////////////////////////////////////

bool TabulateFDistanceOutput(string& strInputFilePathName, string& strAccession, int nOutputCount, string& strOutputFileNameSuffix, string& strTableEntry, string& strErrorEntry)
{
	//  Base file path name
	string strBaseFilePathName = "";
	//  Individual F-Distance file path name
	string strFDFilePathName = "";
	//  Individual file text
	string strFDFileText = "";
	//  F-Distance mutation strings
	vector<string> vMutationCounts;
	//  Path delimiter
	char chrPathDelimiter = strPathDelimiter;

	try
	{
		//  If output table file path name list is not empty
		if (!strInputFilePathName.empty())
		{
			if (!strAccession.empty())
			{
				//  Get base file path name
				strBaseFilePathName = GetBasePath(strInputFilePathName) + chrPathDelimiter + GetFileNameExceptLastExtension(strInputFilePathName);

				//  Iterate output files and tabulate scores
				for (int nCount = -1; nCount < nOutputCount; nCount++)
				{
					if (nCount < 0)
						strFDFilePathName = strBaseFilePathName + "." + strOutputFileNameSuffix + ".fdist";
					else
						strFDFilePathName = strBaseFilePathName + ".sh_" + ConvertIntToString(nCount) + "." + strOutputFileNameSuffix + ".fdist";

					//  Get F-Distance file text
					if (GetFileText(strFDFilePathName, strFDFileText))
					{
						//  If F-Distance file text is not empty
						if (!strFDFileText.empty())
						{
							//  Get F-Distance output mutation count strings
							SplitStringAllowEmptyEntries(strFDFileText, "\n", vMutationCounts);

							//  If mutation count strings are not empty
							if (vMutationCounts.size() > 0)
							{
								if (nCount < 0)
								{
									if (vMutationCounts.size() == 2)
										CompileFDistanceTableOutput(strAccession, vMutationCounts[0], vMutationCounts[1], strTableEntry, true, false);
									else
									{
										string strEmpty = "";

										CompileFDistanceTableOutput(strAccession, vMutationCounts[0], strEmpty, strTableEntry, true, false);
									}
								}
								else
								{
									if (vMutationCounts.size() == 2)
										CompileFDistanceTableOutput(strAccession, vMutationCounts[0], vMutationCounts[1], strTableEntry, true, true);
									else
									{
										string strEmpty = "";

										CompileFDistanceTableOutput(strAccession, vMutationCounts[0], strEmpty, strTableEntry, true, true);
									}
								}
							}
							else
							{
								strErrorEntry = strFDFilePathName + "~Empty Muation Count String(s) [" + ConvertIntToString(nCount) + "]\n";

								ReportTimeStamp("[TabulateFDistanceOutput]", "ERROR:  Input File [" + strFDFilePathName + "] [" + ConvertIntToString(nCount) + "] Mutation Count Container is Empty");
							}
						}
						else
						{
							strErrorEntry = strFDFilePathName + "~Empty File Text\n";

							ReportTimeStamp("[TabulateFDistanceOutput]", "ERROR:  Input File [" + strFDFilePathName + "] [" + ConvertIntToString(nCount) + "] is Empty");
						}
					}
					else
					{
						strErrorEntry = strInputFilePathName + "~Open Failed\n";

						ReportTimeStamp("[TabulateFDistanceOutput]", "ERROR:  Input File [" + strFDFilePathName + "] [" + ConvertIntToString(nCount) + "] Open Failed");
					}
				}

				return true;
			}
			else
			{
				strErrorEntry = strInputFilePathName + "~Empty Accession\n";

				ReportTimeStamp("[TabulateFDistanceOutput]", "ERROR:  Input File [" + strInputFilePathName + "] Accession is Empty");
			}
		}
		else
		{
			ReportTimeStamp("[TabulateFDistanceOutput]", "ERROR:  Input File Path Name is Empty");
		}
	}
	catch (exception ex)
	{
		cout << "ERROR [TabulateFDistanceOutput] Exception Code:  " << ex.what() << "\n";
	}

	return false;
}

////////////////////////////////////////////////////////////////////////////////
//
//  Contabulates BIG F-Distance output from genomic and shuffled F-Distance output files, combining them with the foreground genome's catalog
//         
////////////////////////////////////////////////////////////////////////////////
//
//  [string] strInputFilePathNameList     :  input file path name list
//  [string] strInputFilePathNameTransform:  input file path name transform (includes string replacements, see help)
//  [string] strOutputTableFilePathName   :  output table file path name
//  [string] strCatalogFilePathName       :  foreground catalog file path name
//  [long] lMaxCatalogSize                :  maximum foreground catalog size
//  [int] nOutputCount                    :  original Shuffler output file count
//  [string&] strOutputFileNameSuffix     :  original F-Distance output file name suffix
//  [string] strErrorFilePathName         :  error file path name
//  [int] nMaxProcs                       :  maximum processors for openMP
//                                        :  returns true, if successful; else, false
//         
////////////////////////////////////////////////////////////////////////////////

bool ListContabulateFDistanceOutput(string strInputFilePathNameList, string strInputFilePathNameTransform, string strOutputTableFilePathName, string strCatalogFilePathName, long lMaxCatalogSize, int nOutputCount,
	string strOutputFileNameSuffix, string strErrorFilePathName, int nMaxProcs)
{
	//  List lock
	omp_lock_t lockList;
	//  List file text
	string strInputListFileText = "";
	//  File path names extracted from the file path name list file
	vector<string> vFilePathNames;
	//  Output table entries
	vector<string> vOutputTableEntries;
	//  Output table file text
	string strOutputTableFileText = "";
	//  Error file text entries
	vector<string> vErrorEntries;
	//  Error file text
	string strErrorFileText = "";

	try
	{
		//  If output table file path name list is not empty
		if (!strOutputTableFilePathName.empty())
		{
			//  If input file path name list is not empty
			if (!strInputFilePathNameList.empty())
			{
				//  If catalog file path name list is not empty
				if (!strCatalogFilePathName.empty())
				{
					//  If catalog size is greater than zero
					if (lMaxCatalogSize > 0)
					{
						//  Catalog
						CBase123_Catalog bCatalog(lMaxCatalogSize);

						//  Open catalog
						if (bCatalog.OpenCatalog(strCatalogFilePathName))
						{
							//  Get list file text
							if (GetFileText(strInputFilePathNameList, strInputListFileText))
							{
								//  Split file text into file list
								SplitString(strInputListFileText, '\n', vFilePathNames);

								//  If file path names exist
								if (vFilePathNames.size() > 0)
								{
									//  Initialize output table entries
									vOutputTableEntries.resize(vFilePathNames.size());

									//  Initialize error file entries
									vErrorEntries.resize(vFilePathNames.size());

									//  Get table header
									strOutputTableFileText = bCatalog.GetDemographicsHeader() + GetContabulatedFDistanceOutputTableHeader(nOutputCount);

									//  Initialize time stamp lock
									omp_init_lock(&lockList);

									//  Declare team size
									#pragma omp parallel num_threads(nMaxProcs)
									{
										//  Iterate and process files
										#pragma omp for
										for (long lCount = 0; lCount < vFilePathNames.size(); lCount++)
										{
											//  Test max procs
											if (lCount == 0)
											{
												omp_set_lock(&lockList);
												ReportTimeStamp("[ListContabulateFDistanceOutput]", "NOTE:  Thread Count = " + ConvertIntToString(omp_get_num_threads()));
												omp_unset_lock(&lockList);
											}

											//  Update for timestamp every 10,000 files
											if (lCount % 10000 == 0)
											{
												omp_set_lock(&lockList);
												ReportTimeStamp("[ListContabulateFDistanceOutput]", "NOTE:  Processing Entry [" + ConvertLongToString(lCount) + "] [" + vFilePathNames[lCount] + "]");
												omp_unset_lock(&lockList);
											}

											//  If the file path name is not empty, then proceed
											if (!vFilePathNames[lCount].empty())
											{
												//  Working file path name
												string strWorkingFilePathName = "";

												//  If input file path name transform is not empty
												if (!strInputFilePathNameTransform.empty())
													strWorkingFilePathName = TransformFilePathName(vFilePathNames[lCount], strInputFilePathNameTransform, "");
												else
													strWorkingFilePathName = vFilePathNames[lCount];

												//  If the file is in BIG .fa format
												if (strWorkingFilePathName.substr(strWorkingFilePathName.length() - 3, 3) == ".fa")
												{
													//  Accession
													string strAccession = "";

													//  Get accession
													strAccession = GetAccessionFromBIGFilePathName(strWorkingFilePathName);

													//  If accession is not empty
													if (!strAccession.empty())
													{
														CBase123_Catalog_Entry bceGet;

														//  If catalog search succeeds
														if (bCatalog.GetEntryByAccession(strAccession, bceGet))
														{
															//  Get entry demographics
															vOutputTableEntries[lCount] = bceGet.GetDemographics();

															//  If entry demographics is not empty
															if (!vOutputTableEntries[lCount].empty())
															{
																//  Tabulate this file
																if (!TabulateFDistanceOutput(strWorkingFilePathName, strAccession, nOutputCount, strOutputFileNameSuffix, vOutputTableEntries[lCount], vErrorEntries[lCount]))
																{
																	vErrorEntries[lCount] = strWorkingFilePathName + "~Tabulation Failed\n";

																	omp_set_lock(&lockList);
																	ReportTimeStamp("[ListContabulateFDistanceOutput]", "ERROR:  [" + strWorkingFilePathName + "] Tabulation Failed");
																	omp_unset_lock(&lockList);
																}
															}
															else
															{
																vErrorEntries[lCount] = strWorkingFilePathName + "~Empty Catalog Demographics\n";

																omp_set_lock(&lockList);
																ReportTimeStamp("[ListContabulateFDistanceOutput]", "ERROR:  [" + strWorkingFilePathName + "] Catalog Demographics Entry is Empty");
																omp_unset_lock(&lockList);
															}
														}
														else
														{
															vErrorEntries[lCount] = strWorkingFilePathName + "~Catalog Search Error\n";

															omp_set_lock(&lockList);
															ReportTimeStamp("[ListContabulateFDistanceOutput]", "ERROR:  [" + strWorkingFilePathName + "] Catalog Search Failed");
															omp_unset_lock(&lockList);
														}
													}
													else
													{
														vErrorEntries[lCount] = strWorkingFilePathName + "~Empty Accession\n";

														omp_set_lock(&lockList);
														ReportTimeStamp("[ListContabulateFDistanceOutput]", "ERROR:  [" + strWorkingFilePathName + "] Accession is Empty");
														omp_unset_lock(&lockList);
													}
												}
												else
												{
													vErrorEntries[lCount] = strWorkingFilePathName + "~Format Error\n";

													omp_set_lock(&lockList);
													ReportTimeStamp("[ListContabulateFDistanceOutput]", "ERROR:  Input File [" + strWorkingFilePathName + "] is Not in BIG .fa Format");
													omp_unset_lock(&lockList);
												}
											}
										}
									}

									//  Destroy time stamp lock
									omp_destroy_lock(&lockList);

									//  Concatenate error and table file text
									for (long lCount = 0; lCount < vFilePathNames.size(); lCount++)
									{
										if (!vOutputTableEntries[lCount].empty())
											strOutputTableFileText += vOutputTableEntries[lCount] + "\n";

										if (!vErrorEntries[lCount].empty())
											strErrorFileText += vErrorEntries[lCount] + "\n";
									}

									//  If output table file path name is not empty
									if (!strOutputTableFilePathName.empty())
									{
										//  Write the file
										if (!WriteFileText(strOutputTableFilePathName, strOutputTableFileText))
											ReportTimeStamp("[ListContabulateFDistanceOutput]", "ERROR:  Output Table [" + strOutputTableFilePathName + "] Write Failed");
									}

									//  If error file path name and text are not empty
									if (!strErrorFilePathName.empty())
									{
										//  Write the file
										if (!WriteFileText(strErrorFilePathName, strErrorFileText))
											ReportTimeStamp("[ListContabulateFDistanceOutput]", "ERROR:  Error File [" + strErrorFilePathName + "] Write Failed");
									}

									return true;
								}
								else
								{
									ReportTimeStamp("[ListContabulateFDistanceOutput]", "ERROR:  Input List File [" + strInputFilePathNameList + "] File Text is Empty");
								}
							}
							else
							{
								ReportTimeStamp("[ListContabulateFDistanceOutput]", "ERROR:  Input List File [" + strInputFilePathNameList + "] Open Failed");
							}
						}
						else
						{
							ReportTimeStamp("[ListContabulateFDistanceOutput]", "ERROR:  Catalog [" + strCatalogFilePathName + "] Open Failed");
						}
					}
					else
					{
						ReportTimeStamp("[ListContabulateFDistanceOutput]", "ERROR:  Maximum Catalog Size Must be Greater Than Zero");
					}
				}
				else
				{
					ReportTimeStamp("[ListContabulateFDistanceOutput]", "ERROR:  Catalog File Path Name is Empty");
				}
			}
			else
			{
				ReportTimeStamp("[ListContabulateFDistanceOutput]", "ERROR:  Input List File Path Name is Empty");
			}
		}
		else
		{
			ReportTimeStamp("[ListContabulateFDistanceOutput]", "ERROR:  Output Table File Path Name is Empty");
		}
	}
	catch (exception ex)
	{
		cout << "ERROR [ListContabulateFDistanceOutput] Exception Code:  " << ex.what() << "\n";
	}

	return false;
}

////////////////////////////////////////////////////////////////////////////////
//
//  Clears BIG F-Distance output
//         
////////////////////////////////////////////////////////////////////////////////
//
//  [string] strInputFilePathNameList     :  input file path name list
//  [string] strInputFilePathNameTransform:  input file path name transform (includes string replacements, see help)
//  [int] nOutputCount                    :  original Shuffler output file count
//  [string&] strOutputFileNameSuffix     :  original F-Distance output file name suffix
//  [string] strErrorFilePathName         :  error file path name
//  [int] nMaxProcs                       :  maximum processors for openMP
//                                        :  returns true, if successful; else, false
//         
////////////////////////////////////////////////////////////////////////////////

bool ListClearFDistanceOutput(string strInputFilePathNameList, string strInputFilePathNameTransform, int nOutputCount, string strOutputFileNameSuffix, string strErrorFilePathName, int nMaxProcs)
{
	//  List lock
	omp_lock_t lockList;
	//  List file text
	string strInputListFileText = "";
	//  File path names extracted from the file path name list file
	vector<string> vFilePathNames;
	//  Error file text entries
	vector<string> vErrorEntries;
	//  Error file text
	string strErrorFileText = "";

	try
	{
		//  If input file path name list is not empty
		if (!strInputFilePathNameList.empty())
		{
			//  Get list file text
			if (GetFileText(strInputFilePathNameList, strInputListFileText))
			{
				//  Split file text into file list
				SplitString(strInputListFileText, '\n', vFilePathNames);

				//  If file path names exist
				if (vFilePathNames.size() > 0)
				{
					//  Initialize error file entries
					vErrorEntries.resize(vFilePathNames.size());

					//  Initialize time stamp lock
					omp_init_lock(&lockList);

					//  Declare team size
					#pragma omp parallel num_threads(nMaxProcs)
					{
						//  Iterate and process files
						#pragma omp for
						for (long lCount = 0; lCount < vFilePathNames.size(); lCount++)
						{
							//  Test max procs
							if (lCount == 0)
							{
								omp_set_lock(&lockList);
								ReportTimeStamp("[ListClearFDistanceOutput]", "NOTE:  Thread Count = " + ConvertIntToString(omp_get_num_threads()));
								omp_unset_lock(&lockList);
							}

							//  Update for timestamp every 10,000 files
							if (lCount % 10000 == 0)
							{
								omp_set_lock(&lockList);
								ReportTimeStamp("[ListClearFDistanceOutput]", "NOTE:  Processing Entry [" + ConvertLongToString(lCount) + "] [" + vFilePathNames[lCount] + "]");
								omp_unset_lock(&lockList);
							}

							//  Loop through output count to concatenate file name(s)
							for (int nCount = -1; nCount < nOutputCount; nCount++)
							{
								//  Delete file path name
								string strDeleteFilePathName = "";
								//  Working file path name
								string strWorkingFilePathName = "";
								//  Path delimiter
								char chrPathDelimiter = strPathDelimiter;

								//  If input file path name transform is not empty
								if (!strInputFilePathNameTransform.empty())
									strWorkingFilePathName = TransformFilePathName(vFilePathNames[lCount], strInputFilePathNameTransform, "");
								else
									strWorkingFilePathName = vFilePathNames[lCount];

								//  Get base file path name
								strDeleteFilePathName = GetBasePath(strWorkingFilePathName) + chrPathDelimiter + GetFileNameExceptLastExtension(strWorkingFilePathName);

								if (nCount < 0)
								{
									if (strOutputFileNameSuffix.empty())
										strDeleteFilePathName += ".fdist";
									else
										strDeleteFilePathName += "." + strOutputFileNameSuffix + ".fdist";
								}
								else
								{
									if (strOutputFileNameSuffix.empty())
										strDeleteFilePathName += ".sh_" + ConvertIntToString(nCount) + ".fdist";
									else
										strDeleteFilePathName += ".sh_" + ConvertIntToString(nCount) + "." + strOutputFileNameSuffix + ".fdist";
								}

								//  Delete file
								if (!RemoveFile(strDeleteFilePathName))
								{
									//  Update error
									vErrorEntries[lCount] += strDeleteFilePathName + "~Remove Failed\n";

									//  Update console
									omp_set_lock(&lockList);
									ReportTimeStamp("[ListClearFDistanceOutput]", "ERROR:  File [" + strDeleteFilePathName + "] Remove Failed");
									omp_unset_lock(&lockList);
								}
							}
						}
					}

					//  Destroy time stamp lock
					omp_destroy_lock(&lockList);

					//  Concatenate error and table file text
					for (long lCount = 0; lCount < vFilePathNames.size(); lCount++)
					{
						if (!vErrorEntries[lCount].empty())
							strErrorFileText += vErrorEntries[lCount] + "\n";
					}

					//  If error file path name and text are not empty
					if (!strErrorFilePathName.empty())
					{
						//  Write the file
						if (!WriteFileText(strErrorFilePathName, strErrorFileText))
							ReportTimeStamp("[ListClearFDistanceOutput]", "ERROR:  Error File [" + strErrorFilePathName + "] Write Failed");
					}

					return true;
				}
				else
				{
					ReportTimeStamp("[ListClearFDistanceOutput]", "ERROR:  Input List File [" + strInputFilePathNameList + "] File Text is Empty");
				}
			}
			else
			{
				ReportTimeStamp("[ListClearFDistanceOutput]", "ERROR:  Input List File [" + strInputFilePathNameList + "] Open Failed");
			}
		}
		else
		{
			ReportTimeStamp("[ListClearFDistanceOutput]", "ERROR:  Input List File Path Name is Empty");
		}
	}
	catch (exception ex)
	{
		cout << "ERROR [ListClearFDistanceOutput] Exception Code:  " << ex.what() << "\n";
	}

	return false;
}

////////////////////////////////////////////////////////////////////////////////
//
//  Performs the F-Distance analysis on a list of BIG .fa format files
//        
////////////////////////////////////////////////////////////////////////////////
//
//  [string] strOutputTableFilePathName             :  output table file path name
//  [string] strBackgroundFilePathNameList          :  background file path name list (BIG .fa format)
//  [string] strBackgroundInputFilePathNameTransform:  background input file path name transform (includes string replacements, see help)
//  [string] strBackgroundCatalogFilePathName       :  BIG background genome cataloge file path name
//  [long] lMaxBackgroundCatalogSize                :  estimated maximum size of the background catalog
//  [bool] bBackgroundBidirect                      :  process background bidirectional, if true
//  [bool] bBackgroundAllowUnknowns                 :  process background unknown chracters, if true
//  [string] strBackgroundErrorFilePathName         :  base file name of the background error file
//  [string] strForegroundFilePathNameList          :  foreground file path name list (BIG .fa format)
//  [string] strForegroundInputFilePathNameTransform:  foreground input file path name transform (includes string replacements, see help)
//  [string] strForegroundCatalogFilePathName       :  BIG foreground genome cataloge file path name
//  [long] lMaxForegroundCatalogSize                :  estimated maximum size of the foreground catalog
//  [bool] bForegroundBidirect                      :  process foreground bidirectional, if true
//  [bool] bForegroundAllowUnknowns                 :  process foreground unknown chracters, if true
//  [string] strOutputFileNameSuffix                :  output file name suffix
//  [string] strForegroundErrorFilePathName         :  base file name of the foreground error file
//  [int] nNMerLength                               :  nMer length to analyze
//  [int] nMaxProcs                                 :  maximum processor count (for openMP)
//                                                 :  returns true, if successful; else, false
//
////////////////////////////////////////////////////////////////////////////////

bool PerformFDistanceAnalysis(string strOutputTableFilePathName, string strBackgroundFilePathNameList, string strBackgroundInputFilePathNameTransform, string strBackgroundCatalogFilePathName,
	long lMaxBackgroundCatalogSize, bool bBackgroundBidirect, bool bBackgroundAllowUnknowns, string strBackgroundErrorFilePathName, string strForegroundFilePathNameList,
	string strForegroundInputFilePathNameTransform, string strForegroundCatalogFilePathName, long lMaxForegroundCatalogSize, bool bForegroundBidirect, bool bForegroundAllowUnknowns,
	string strOutputFileNameSuffix, string strForegroundErrorFilePathName, int nNMerLength, int nMaxProcs)
{
	//  Return status, is success if true, else is not-error
	bool bStatusSuccess = false;
	//  Output table file text
	vector<string> vOutputTableEntries;

	try
	{
		if (nNMerLength == 8)
			bStatusSuccess = InitializeWriteLock16();
		else if (nNMerLength == 16)
			bStatusSuccess = InitializeWriteLock32();

		//  If write lock initialized
		if(bStatusSuccess)
		{
			//  If output table file path name is not empty
			if (!strOutputTableFilePathName.empty())
			{
				//  If background genome catalog file path name list is not empty
				if (!strBackgroundCatalogFilePathName.empty())
				{
					//  If background genome catalog size is Not zero
					if (lMaxBackgroundCatalogSize > 0)
					{
						//  If background file path name list is not empty
						if (!strBackgroundFilePathNameList.empty())
						{
							//  If foreground genome catalog file path name list is not empty
							if (!strForegroundCatalogFilePathName.empty())
							{
								//  If background genome catalog size is Not zero
								if (lMaxBackgroundCatalogSize > 0)
								{
									//  If foreground file path name list is not empty
									if (!strForegroundFilePathNameList.empty())
									{
										//  If nMer length is properly set
										if (nNMerLength > 0)
										{
											//  If nMaxProcs is properly set
											if (nMaxProcs > 0)
											{
												//  Background genome catalog
												CBase123_Catalog b123BackgroundCatalog(lMaxBackgroundCatalogSize);

												//  Open background catalog
												if (b123BackgroundCatalog.OpenCatalog(strBackgroundCatalogFilePathName))
												{
													//  Update console; end application;
													ReportTimeStamp("[PerformFDistanceAnalysis]", "Background Catalog Opened");

													//  Foreground genome catalog
													CBase123_Catalog b123ForegroundCatalog(lMaxForegroundCatalogSize);

													//  Open foreground catalog
													if (b123ForegroundCatalog.OpenCatalog(strForegroundCatalogFilePathName))
													{
														//  Update console; end application;
														ReportTimeStamp("[PerformFDistanceAnalysis]", "Foreground Catalog Opened");

														//  Initialize the background array
														bStatusSuccess = false;
														if (nNMerLength == 8)
															bStatusSuccess = InitializeBackground16();
														else if (nNMerLength == 16)
															bStatusSuccess = InitializeBackground32();

														//  If background array is set
														if (bStatusSuccess)
														{
															//  Update console; end application;
															ReportTimeStamp("[PerformFDistanceAnalysis]", "Background Collection Initialized");

															//  Destroy the background array
															bStatusSuccess = false;
															if (nNMerLength == 8)
																bStatusSuccess = ProcessFDistanceList16(strBackgroundFilePathNameList, strBackgroundInputFilePathNameTransform, b123BackgroundCatalog, bBackgroundBidirect, nNMerLength, true, bBackgroundAllowUnknowns, strOutputFileNameSuffix, strBackgroundErrorFilePathName, vOutputTableEntries, nMaxProcs);
															else if (nNMerLength == 16)
																bStatusSuccess = ProcessFDistanceList32(strBackgroundFilePathNameList, strBackgroundInputFilePathNameTransform, b123BackgroundCatalog, bBackgroundBidirect, nNMerLength, true, bBackgroundAllowUnknowns, strOutputFileNameSuffix, strBackgroundErrorFilePathName, vOutputTableEntries, nMaxProcs);

															//  Process background file list
															if (bStatusSuccess)
															{
																//  Update console; end application;
																ReportTimeStamp("[PerformFDistanceAnalysis]", "Background Loaded");

																//  Destroy the background array
																bStatusSuccess = false;
																if (nNMerLength == 8)
																	bStatusSuccess = ProcessFDistanceList16(strForegroundFilePathNameList, strForegroundCatalogFilePathName, b123ForegroundCatalog, bForegroundBidirect, nNMerLength, false, bForegroundAllowUnknowns, strOutputFileNameSuffix, strForegroundErrorFilePathName, vOutputTableEntries, nMaxProcs);
																else if (nNMerLength == 16)
																	bStatusSuccess = ProcessFDistanceList32(strForegroundFilePathNameList, strForegroundCatalogFilePathName, b123ForegroundCatalog, bForegroundBidirect, nNMerLength, false, bForegroundAllowUnknowns, strOutputFileNameSuffix, strForegroundErrorFilePathName, vOutputTableEntries, nMaxProcs);

																//  Process foreground file list
																if (bStatusSuccess)
																{
																	//  Update console; end application;
																	ReportTimeStamp("[PerformFDistanceAnalysis]", "Foreground Analyzed");

																	//  Write output table
																	if (!WriteFDistanceOutputTable(strOutputTableFilePathName, vOutputTableEntries))
																	{
																		ReportTimeStamp("[PerformFDistanceAnalysis]", "F-Distance Output Table File [" + strOutputTableFilePathName + "] Write Failed");
																	}
																}
																else
																{
																	ReportTimeStamp("[PerformFDistanceAnalysis]", "ERROR:  Foreground Process Failed");
																}
															}
															else
															{
																ReportTimeStamp("[PerformFDistanceAnalysis]", "ERROR:  Background Background Process Failed");
															}

															//  Destroy the background array
															bStatusSuccess = false;
															if (nNMerLength == 8)
																bStatusSuccess = DestroyBackground16();
															else if (nNMerLength == 16)
																bStatusSuccess = DestroyBackground32();

															//  If error, report
															if (!bStatusSuccess)
																ReportTimeStamp("[PerformFDistanceAnalysis]", "ERROR:  Background Destruction Failed");

															//  Destroy the write-lock
															bStatusSuccess = false;
															if (nNMerLength == 8)
																bStatusSuccess = DestroyWriteLock16();
															else if (nNMerLength == 16)
																bStatusSuccess = DestroyWriteLock32();

															//  If error, report
															if (!bStatusSuccess)
																ReportTimeStamp("[PerformFDistanceAnalysis]", "ERROR:  Write Lock Destruction Failed");

															//  Return success
															return true;
														}
														else
														{
															ReportTimeStamp("[PerformFDistanceAnalysis]", "ERROR:  Background Initialization Failed");
														}

														//  Clear foreground catalog entries
														b123ForegroundCatalog.CloseCatalog();
													}
													else
													{
														ReportTimeStamp("[PerformFDistanceAnalysis]", "ERROR:  Foreground Genome Catalog [" + strForegroundCatalogFilePathName + "] Open Failed");
													}

													//  Clear background catalog entries
													b123BackgroundCatalog.CloseCatalog();
												}
												else
												{
													ReportTimeStamp("[PerformFDistanceAnalysis]", "ERROR:  Background Genome Catalog [" + strBackgroundCatalogFilePathName + "] Open Failed");
												}
											}
											else
											{
												ReportTimeStamp("[PerformFDistanceAnalysis]", "ERROR:  openMP Maximum Processor Count is Not Properly Set:  Should be Greater Than 0");
											}
										}
										else
										{
											ReportTimeStamp("[PerformFDistanceAnalysis]", "ERROR:  nMer Length is Not Properly Set:  Should be 8 or 16");
										}
									}
									else
									{
										ReportTimeStamp("[PerformFDistanceAnalysis]", "ERROR:  Foreground Genome File Path Name List is Empty");
									}
								}
								else
								{
									ReportTimeStamp("[PerformFDistanceAnalysis]", "ERROR:  Foreground Catalog Estimated Maximum Size Must be Greater Than Zero");
								}
							}
							else
							{
								ReportTimeStamp("[PerformFDistanceAnalysis]", "ERROR:  Foreground Genome Catalog File Path Name List is Empty");
							}
						}
						else
						{
							ReportTimeStamp("[PerformFDistanceAnalysis]", "ERROR:  Background Genome File Path Name List is Empty");
						}
					}
					else
					{
						ReportTimeStamp("[PerformFDistanceAnalysis]", "ERROR:  Background Catalog Estimated Maximum Size Must be Greater Than Zero");
					}
				}
				else
				{
					ReportTimeStamp("[PerformFDistanceAnalysis]", "ERROR:  Background Genome Catalog File Path Name List is Empty");
				}
			}
			else
			{
				ReportTimeStamp("[PerformFDistanceAnalysis]", "ERROR:  F-Distance Output Table File Path Name is Empty");
			}
		}
		else
		{
			ReportTimeStamp("[PerformFDistanceAnalysis]", "ERROR:  Write-Lock Initialization Failed");
		}
	}
	catch (exception ex)
	{
		cout << "ERROR [PerformFDistanceAnalysis] Exception Code:  " << ex.what() << "\n";
	}

	return false;
}
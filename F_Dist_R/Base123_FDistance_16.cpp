// Base123_FDistance_16.cpp : Performs Base123 Fofanov Distance Genomic Analysis

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
//  6 July 2016
//  (All Rights Reserved)
//
////////////////////////////////////////////////////////////////////////////////

#include "F_Dist_R.h"
#include "Base123_Catalog_Entry.h"
#include "Base123_Catalog.h"
#include "Base123_FDistance_16.h"
#include "Base123_Utilities.h"

#include <math.h>
#include <sstream>
#include <omp.h>

//  Background array
uint16_t* m_unaBackground16 = NULL;
bool m_bBackgroundPolyTU16 = false;

//  write lock
omp_lock_t writelock16;

////////////////////////////////////////////////////////////////////////////////
//
//  Converts a uint16_t sequence to a string sequence
//         
////////////////////////////////////////////////////////////////////////////////
//
//  [string] strSequence:  sequence to convert and return
//                     :  returns sequence string, if successful; else, an empty string
//         
////////////////////////////////////////////////////////////////////////////////

uint16_t ConvertStringToSequence16(string strSequence)
{
	//  Sequence to return
	uint16_t untSequence = 0;

	try
	{
		//  If background array is set
		if (!strSequence.empty())
		{
			if (strSequence != "uuuuuuuu")
			{
				for (int nCount = 0; nCount < strSequence.length(); nCount++)
				{
					//  Base at this position
					if (strSequence.substr(nCount, 1) == "a")
						untSequence += m_untA;
					else if (strSequence.substr(nCount, 1) == "c")
						untSequence += m_untC;
					else if (strSequence.substr(nCount, 1) == "g")
						untSequence += m_untG;
					else if (strSequence.substr(nCount, 1) == "u")
						untSequence += m_untTU;

					//  Next base
					if (nCount < strSequence.length() - 1)
						untSequence = untSequence << 2;
				}
			}
			else
				untSequence = UINT16_MAX;
		}
		else
		{
			ReportTimeStamp("[ConvertStringToSequence16]", "ERROR:  Sequence is Empty");
		}
	}
	catch (exception ex)
	{
		cout << "ERROR [ConvertStringToSequence16] Exception Code:  " << ex.what() << "\n";
	}

	return untSequence;
}

////////////////////////////////////////////////////////////////////////////////
//
//  Converts a uint16_t sequence to a string sequence
//         
////////////////////////////////////////////////////////////////////////////////
//
//  [uint16_t] untSequence:  sequence to convert and return
//  [int] nNMerLength     :  nMer length to convert
//                       :  returns sequence string, if successful; else, an empty string
//         
////////////////////////////////////////////////////////////////////////////////

string ConvertSequenceToString16(uint16_t untSequence, int nNMerLength)
{
	//  Sequence to return
	string strSequence = "";
	//  Set position mask to last base position
	uint16_t untPositionMask = 0b1100000000000000;

	try
	{
		//  If background array is set
		if ((untSequence >= 0) && (untSequence < UINT16_MAX))
		{
			for (int nCount = nNMerLength - 1; nCount >= 0; nCount--)
			{
				//  Base at this position
				uint16_t untThisBase = untSequence & untPositionMask;
				untThisBase = untThisBase >> nCount * 2;

				if (untThisBase == m_untA)
					strSequence += "a";
				else if (untThisBase == m_untC)
					strSequence += "c";
				else if (untThisBase == m_untG)
					strSequence += "g";
				else if (untThisBase == m_untTU)
					strSequence += "u";

				//  Next base
				untPositionMask = untPositionMask >> 2;
			}
		}
		else if (untSequence == UINT16_MAX)
			strSequence = "uuuuuuuu";
		else
		{
			ReportTimeStamp("[ConvertSequenceToString16]", "ERROR:  Sequence [" + ConvertUnsignedInt16ToString(untSequence) + "] is Out of Range [0:" + ConvertUnsignedInt16ToString(UINT16_MAX) + "]");
		}
	}
	catch (exception ex)
	{
		cout << "ERROR [ConvertSequenceToString16] Exception Code:  " << ex.what() << "\n";
	}

	return strSequence;
}

////////////////////////////////////////////////////////////////////////////////
//
//  Writes the background array to file
//         
////////////////////////////////////////////////////////////////////////////////
//
//  [string] strOutputFilePathBaseName:  output file path name
//  [int] nNMerLength                 :  nMer length to convert
//                                   :  returns true, if successful; else, false
//         
////////////////////////////////////////////////////////////////////////////////

bool WriteBackgroundArray16(string strOutputFilePathBaseName, int nNMerLength)
{
	//  Present file test
	string strPresentFileText = "";
	//  Absent file text
	string strAbsentFileText = "";

	try
	{
		//  If background array is set
		if (m_unaBackground16 != NULL)
		{
			//  If output file path base name is not empty
			if (!strOutputFilePathBaseName.empty())
			{
				//  Iterate background and concatenate file texts
				for (uint16_t untCount = 0; untCount < UINT16_MAX; untCount++)
				{
					if (m_unaBackground16[untCount] == 0)
						strAbsentFileText += ConvertSequenceToString16(untCount, nNMerLength) + "\n";
					else
						strPresentFileText += ConvertSequenceToString16(untCount, nNMerLength) + "\n";
				}

				//  Write output
				if (!WriteFileText(strOutputFilePathBaseName + "_absent.txt", strAbsentFileText))
					ReportTimeStamp("[WriteBackgroundArray16]", "ERROR:  Absent File Write Failed");

				if (!WriteFileText(strOutputFilePathBaseName + "_present.txt", strPresentFileText))
					ReportTimeStamp("[WriteBackgroundArray16]", "ERROR:  Present File Write Failed");
			}
			else
			{
				ReportTimeStamp("[WriteBackgroundArray16]", "ERROR:  Output File Path Base Name is Empty");
			}
		}
		else
		{
			ReportTimeStamp("[WriteBackgroundArray16]", "ERROR:  Background Collection is Not Set");
		}
	}
	catch (exception ex)
	{
		cout << "ERROR [WriteBackgroundArray16] Exception Code:  " << ex.what() << "\n";
	}

	return false;
}

////////////////////////////////////////////////////////////////////////////////
//
//  Mutates the sequence one base at a time
//         
////////////////////////////////////////////////////////////////////////////////
//
//  [uint16_t&] untSequence:  sequence to mutate
//                        :  returns true, if sequence is found in background with one mutation; else, false
//         
////////////////////////////////////////////////////////////////////////////////

bool MutateOne16(uint16_t& untSequence)
{
	try
	{
		//  Set position untPositionMask1 to last base position
		uint16_t untPositionMask1 = 0b11;

		//  If background array is set
		if (m_unaBackground16 != NULL)
		{
			//  Iterate through base positions
			while (untPositionMask1 != 0)
			{
				//  Base at this position to be mutated
				uint16_t untThisBase = untSequence & untPositionMask1;

				//  Mutation base
				uint16_t untMutationBase = untThisBase;

				//  Sequence untPositionMask1 with position hole
				uint16_t untSequenceMask = 0b1111111111111111 & untPositionMask1;

				//  Iterate every other base at this position
				//  (series at += is:  00 & 11 -> 11 & 11 -> 10 & 11 -> 01 & 11 -> 00 (original base)
				//  So anding 11 three times gets ever possible next base EXCEPT current base
				for (int lCountBases = 0; lCountBases < 3; lCountBases++)
				{
					//  Add mutation base to sequence untPositionMask1 to get next base at this position
					untMutationBase += untSequenceMask;

					//  Mask out this next base
					untMutationBase = untMutationBase & untPositionMask1;

					//  Flip untPositionMask1 bits to make a positional hole for insertaton of this next base
					uint16_t untInsertMask = ~untPositionMask1;

					//  Insert the new base at the current position in the original sequence to set the next mutational index
					uint16_t untIndex = (untSequence & untInsertMask) + untMutationBase;

					//  Sequence is found in background after one mutation, any position
					if (untSequence == UINT16_MAX)
					{
						if (m_bBackgroundPolyTU16)
							return true;
					}
					else if ((untIndex >= 0) && (untIndex < UINT16_MAX) && (m_unaBackground16[untIndex] == 1))
						return true;
				}

				//  Shift position untPositionMask1 leftward to next base position
				untPositionMask1 = untPositionMask1 << 2;
			}
		}
		else
		{
			ReportTimeStamp("[MutateOne16]", "ERROR:  Background Collection is Not Set");
		}
	}
	catch (exception ex)
	{
		cout << "ERROR [MutateOne16] Exception Code:  " << ex.what() << "\n";
	}

	return false;
}

////////////////////////////////////////////////////////////////////////////////
//
//  Mutates the sequence two base(s) at a time
//         
////////////////////////////////////////////////////////////////////////////////
//
//  [uint16_t&] untSequence:  sequence to mutate
//                        :  returns true, if sequence is found in background with one mutation; else, false
//         
////////////////////////////////////////////////////////////////////////////////

bool MutateTwo16(uint16_t& untSequence)
{
	try
	{
		//  Set position untPositionMask1 to last base position
		uint16_t untPositionMask1 = 0b11;
		//  Mutation base 1
		uint16_t untMutationBase1 = untPositionMask1;
		//  Set second position untPositionMask1 to net to last base position
		uint16_t untPositionMask2 = 0b1100;
		//  Mutation base 2
		uint16_t untMutationBase2 = untPositionMask2;
		//  Combined untPositionMask1
		uint16_t untCombinedMask = 0;

		//  If background array is set
		if (m_unaBackground16 != NULL)
		{
			//  Iterate through base positions
			while (untPositionMask1 != 0)
			{
				//  Iterate every base at second position
				for (int lCountBases2 = 0; lCountBases2 < 4; lCountBases2++)
				{
					//  Get next base at position 2
					untMutationBase2 += untPositionMask2;
					//  Isolate next base at position 2
					untMutationBase2 = untMutationBase2 & untPositionMask2;

					//  Iterate every base at first position
					for (int lCountBases1 = 0; lCountBases1 < 4; lCountBases1++)
					{
						//  Get next base at position 1
						untMutationBase1 += untPositionMask1;
						//  Isolate next base at position 1
						untMutationBase1 = untMutationBase1 & untPositionMask1;

						//  Combine the two mutation bases
						uint16_t untCombinedMutationMask = untMutationBase2 + untMutationBase1;

						//  Cbomine the position untPositionMask1s
						untCombinedMask = untPositionMask2 + untPositionMask1;

						//  Concatanet mutated sequence
						uint16_t untIndex = (untSequence & ~untCombinedMask) + untCombinedMutationMask;

						//  Sequence is found in background after two mutations, any positions
						if (untSequence == UINT16_MAX)
						{
							if (m_bBackgroundPolyTU16)
								return true;
						}
						else if ((untIndex >= 0) && (untIndex < UINT16_MAX) && (m_unaBackground16[untIndex] == 1))
							return true;
					}
				}

				//  Move to next mutational base position 2
				untPositionMask2 = untPositionMask2 << 2;
				//  If mutational base position 2 is at beginning of sequence, reset to the end
				if (untPositionMask2 == 0)
					untPositionMask2 = 0b11;

				//  If mutational base position 2 == mutational base position 1, then shift both left to start next position 1
				if (untPositionMask2 == untPositionMask1)
				{
					untPositionMask1 = untPositionMask1 << 2;
					untPositionMask2 = untPositionMask1 << 2;
				}
			}
		}
		else
		{
			ReportTimeStamp("[MutateTwo16]", "ERROR:  Background Collection is Not Set");
		}
	}
	catch (exception ex)
	{
		cout << "ERROR [MutateTwo16] Exception Code:  " << ex.what() << "\n";
	}

	return false;
}

////////////////////////////////////////////////////////////////////////////////
//
//  Mutates the sequence three base(s) at a time
//         
////////////////////////////////////////////////////////////////////////////////
//
//  [uint16_t&] untSequence:  sequence to mutate
//                        :  returns true, if sequence is found in background with one mutation; else, false
//         
////////////////////////////////////////////////////////////////////////////////

bool MutateThree16(uint16_t& untSequence)
{
	try
	{
		uint16_t untPositionMask1 = 0b11;
		uint16_t untMutationBase1 = untPositionMask1;
		uint16_t untPositionMask2 = 0b1100;
		uint16_t untMutationBase2 = untPositionMask2;
		uint16_t untPositionMask3 = 0b110000;
		uint16_t untMutationBase3 = untPositionMask3;
		uint16_t untCombinedMask = 0;

		//  If background array is set
		if (m_unaBackground16 != NULL)
		{
			while (untPositionMask1 != 0)
			{
				for (int nCount3 = 0; nCount3 < 4; nCount3++)
				{
					untMutationBase3 += untPositionMask3;
					untMutationBase3 = untMutationBase3 & untPositionMask3;
					for (int nCount2 = 0; nCount2 < 4; nCount2++)
					{
						untMutationBase2 += untPositionMask2;
						untMutationBase2 = untMutationBase2 & untPositionMask2;
						for (int nCount1 = 0; nCount1 < 4; nCount1++)
						{
							untMutationBase1 += untPositionMask1;
							untMutationBase1 = untMutationBase1 & untPositionMask1;
							uint16_t untCombinedMutationMask = untMutationBase3 + untMutationBase2 + untMutationBase1;
							untCombinedMask = untPositionMask3 + untPositionMask2 + untPositionMask1;
							uint16_t untIndex = (untSequence & ~untCombinedMask) + untCombinedMutationMask;
							if (untSequence == UINT16_MAX)
							{
								if (m_bBackgroundPolyTU16)
									return true;
							}
							else if ((untIndex >= 0) && (untIndex < UINT16_MAX) && (m_unaBackground16[untIndex] == 1))
								return true;
						}
					}
				}
				untPositionMask3 = untPositionMask3 << 2;
				if (untPositionMask3 == 0)
				{
					untPositionMask2 = untPositionMask2 << 2;
					if (untPositionMask2 == 0)
					{
						untPositionMask1 = untPositionMask1 << 2;
						untPositionMask2 = untPositionMask1 << 2;
					}
					untPositionMask3 = untPositionMask2 << 2;
				}
				untCombinedMask = untPositionMask1 | untPositionMask2 | untPositionMask3;
			}
		}
		else
		{
			ReportTimeStamp("[MutateThree16]", "ERROR:  Background Collection is Not Set");
		}
	}
	catch (exception ex)
	{
		cout << "ERROR [MutateThree16] Exception Code:  " << ex.what() << "\n";
	}

	return false;
}

////////////////////////////////////////////////////////////////////////////////
//
//  Mutates the sequence four base(s) at a time
//         
////////////////////////////////////////////////////////////////////////////////
//
//  [uint16_t&] untSequence:  sequence to mutate
//                        :  returns true, if sequence is found in background with one mutation; else, false
//         
////////////////////////////////////////////////////////////////////////////////

bool MutateFour16(uint16_t& untSequence)
{
	try
	{
		uint16_t untPositionMask1 = 0b11;
		uint16_t untMutationBase1 = untPositionMask1;
		uint16_t untPositionMask2 = 0b1100;
		uint16_t untMutationBase2 = untPositionMask2;
		uint16_t untPositionMask3 = 0b110000;
		uint16_t untMutationBase3 = untPositionMask3;
		uint16_t untPositionMask4 = 0b11000000;
		uint16_t untMutationBase4 = untPositionMask4;
		uint16_t untCombinedMask = 0;

		//  If background array is set
		if (m_unaBackground16 != NULL)
		{
			while (untPositionMask1 != 0)
			{
				for (int nCount4 = 0; nCount4 < 4; nCount4++)
				{
					untMutationBase4 += untPositionMask4;
					untMutationBase4 = untMutationBase4 & untPositionMask4;
					for (int nCount3 = 0; nCount3 < 4; nCount3++)
					{
						untMutationBase3 += untPositionMask3;
						untMutationBase3 = untMutationBase3 & untPositionMask3;
						for (int nCount2 = 0; nCount2 < 4; nCount2++)
						{
							untMutationBase2 += untPositionMask2;
							untMutationBase2 = untMutationBase2 & untPositionMask2;
							for (int nCount1 = 0; nCount1 < 4; nCount1++)
							{
								untMutationBase1 += untPositionMask1;
								untMutationBase1 = untMutationBase1 & untPositionMask1;
								uint16_t untCombinedMutationMask = untMutationBase4 + untMutationBase3 + untMutationBase2 + untMutationBase1;
								untCombinedMask = untPositionMask4 + untPositionMask3 + untPositionMask2 + untPositionMask1;
								uint16_t untIndex = (untSequence & ~untCombinedMask) + untCombinedMutationMask;
								if (untSequence == UINT16_MAX)
								{
									if (m_bBackgroundPolyTU16)
										return true;
								}
								else if ((untIndex >= 0) && (untIndex < UINT16_MAX) && (m_unaBackground16[untIndex] == 1))
									return true;
							}
						}
					}
				}
				untPositionMask4 = untPositionMask4 << 2;
				if (untPositionMask4 == 0)
				{
					untPositionMask3 = untPositionMask3 << 2;
					if (untPositionMask3 == 0)
					{
						untPositionMask2 = untPositionMask2 << 2;
						if (untPositionMask2 == 0)
						{
							untPositionMask1 = untPositionMask1 << 2;
							untPositionMask2 = untPositionMask1 << 2;
						}
						untPositionMask3 = untPositionMask2 << 2;
					}
					untPositionMask4 = untPositionMask3 << 2;
				}
			}
		}
		else
		{
			ReportTimeStamp("[MutateFour16]", "ERROR:  Background Collection is Not Set");
		}
	}
	catch (exception ex)
	{
		cout << "ERROR [MutateFour16] Exception Code:  " << ex.what() << "\n";
	}

	return false;
}

////////////////////////////////////////////////////////////////////////////////
//
//  Mutates the sequence until it is found in the background or until mutational limit is exceeded
//         
////////////////////////////////////////////////////////////////////////////////
//
//  [uint16_t&] untSequence:  sequence to mutate
//                        :  returns mutational count, if successful; else, -1
//         
////////////////////////////////////////////////////////////////////////////////

int MutateSequence16(uint16_t& untSequence)
{
	try
	{
		//  If background array is set
		if (m_unaBackground16 != NULL)
		{
			//  Get mutational count and store it at first position
			if ((untSequence == UINT16_MAX) && (m_bBackgroundPolyTU16))
				return 0;

			if (m_unaBackground16[untSequence] == 1)
				return 0;

			if (MutateOne16(untSequence))
				return 1;

			if (MutateOne16(untSequence))
				return 1;

			if (MutateTwo16(untSequence))
				return 2;

			if (MutateThree16(untSequence))
				return 3;

			if (MutateFour16(untSequence))
				return 4;
		}
		else
		{
			ReportTimeStamp("[MutateSequence16]", "ERROR:  Background Collection is Not Set");
		}
	}
	catch (exception ex)
	{
		cout << "ERROR [MutateSequence16] Exception Code:  " << ex.what() << "\n";
	}

	return -1;
}

////////////////////////////////////////////////////////////////////////////////
//
//  Processes the foreground genome; identifies foreground nMers present in the background
//      or performs mutational analysis until nMer is found;
//         
////////////////////////////////////////////////////////////////////////////////
//
//  [string&] strOutputFilePathName  :  output file path name
//  [string&] strAccession           :  accession of sequence to process
//  [string&] strSequence            :  sequence to process
//  [bool] bBidirectional            :  process bidirectionally, if true
//  [bool] bForegroundAllowUnknowns  :  process foreground unknown chracters, if true
//  [int] nNMerLength                :  nMer length to analyze
//  [string&] strOutputTableEntry    :  F-Distance table file text to concatenate
//                                  :  returns true, if successful; else, false
//         
////////////////////////////////////////////////////////////////////////////////

bool ProcessForeground16(string& strOutputFilePathName, string& strAccession, string& strSequence, bool bBidirectional, bool bForegroundAllowUnknowns, int nNMerLength, 
	string& strOutputTableEntry)
{
	//  Character sequence
	string strSubSequence = "";
	//  nMer sequence, binary
	uint16_t untSubSequence = 0;
	//  Mutational count
	int nMutationCount = 0;
	//  Process iteration maximum, according to bidirectional flag
	int nMaxProcess = 1;
	//  Forward output string
	string strForwardOutput = "";
	//  Reverse output string
	string strReverseOutput = "";
	//  Output file text
	string strOutputFileText = "";

	try
	{
		//  Foreground nMer to mutate
		uint16_t untNMer = 0b00;

		//  If background array is set
		if (m_unaBackground16 != NULL)
		{
			//  If accession is not empty
			if (!strAccession.empty())
			{
				//  If output file path name is not empty
				if (!strOutputFilePathName.empty())
				{
					//  If the input sequence is not empty
					if (!strSequence.empty())
					{
						//  Reject sequences with gaps
						if (strSequence.find_first_of('-') != string::npos)
						{
							ReportTimeStamp("[ProcessForeground16]", "ERROR:  Foreground Sequence Contains a Gap of Indeterminate Length");

							return false;
						}
						//  Reject sequences with poly-n
						else if (strSequence.find("nnnnnnnn") != string::npos)
						{
							ReportTimeStamp("[ProcessForeground16]", "ERROR:  Input Sequence Contains Excessive Poly-n Bases");

							return false;
						}

						//  If bidirectional, iterate process twice
						if (bBidirectional)
							nMaxProcess = 2;

						//  Process iteration maximum, according to bidirectional flag
						for (int nCountPass = 0; nCountPass < nMaxProcess; nCountPass++)
						{
							//  Get reverse complement on second pass
							if (nCountPass > 0)
								strSequence = ConvertToReverseCompliment(strSequence);

							//  Get first nMer
							strSubSequence = strSequence.substr(0, nNMerLength);

							//  Process known bases
							if ((strSubSequence.find_first_of('r') != string::npos) ||
								(strSubSequence.find_first_of('y') != string::npos) ||
								(strSubSequence.find_first_of('k') != string::npos) ||
								(strSubSequence.find_first_of('m') != string::npos) ||
								(strSubSequence.find_first_of('s') != string::npos) ||
								(strSubSequence.find_first_of('w') != string::npos) ||
								(strSubSequence.find_first_of('b') != string::npos) ||
								(strSubSequence.find_first_of('d') != string::npos) ||
								(strSubSequence.find_first_of('h') != string::npos) ||
								(strSubSequence.find_first_of('v') != string::npos) ||
								(strSubSequence.find_first_of('n') != string::npos))
							{
								if (bForegroundAllowUnknowns)
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

									strSubSequence = ReplaceInString(strSubSequence, "r", "a", false);
									strSubSequence = ReplaceInString(strSubSequence, "y", "c", false);
									strSubSequence = ReplaceInString(strSubSequence, "k", "g", false);
									strSubSequence = ReplaceInString(strSubSequence, "m", "a", false);
									strSubSequence = ReplaceInString(strSubSequence, "s", "c", false);
									strSubSequence = ReplaceInString(strSubSequence, "w", "a", false);
									strSubSequence = ReplaceInString(strSubSequence, "b", "c", false);
									strSubSequence = ReplaceInString(strSubSequence, "d", "a", false);
									strSubSequence = ReplaceInString(strSubSequence, "h", "a", false);
									strSubSequence = ReplaceInString(strSubSequence, "v", "a", false);
									strSubSequence = ReplaceInString(strSubSequence, "n", "a", false);
								}
								else
								{
									ReportTimeStamp("[ProcessForeground16]", "ERROR:  Foreground Sequence Contains Unknown Characters and -background_allow_unknowns [-bau] is Not Set");

									return false;
								}
							}

							//  Get binary sequence
							untSubSequence = ConvertStringToSequence16(strSubSequence);

							//  Get mutation count
							nMutationCount = MutateSequence16(untSubSequence);
							if (nMutationCount >= 0)
							{
								//  Concatenate forward output
								if (nCountPass == 0)
									strForwardOutput += ConvertIntToString(nMutationCount);
								else
									strReverseOutput += ConvertIntToString(nMutationCount);
							}
							else
							{
								ReportTimeStamp("[ProcessForeground16]", "ERROR:  Foreground Sequence [" + ConvertSequenceToString16(untSubSequence, nNMerLength) + "] @ [0] Mutation Failed or Mutation Count Exceeds Limit [8]");

								return false;
							}

							//  Iterate subsequent nNMerLength characters to build remaing nMers
							for (long lCountBases = nNMerLength; lCountBases < (strSequence.length() - nNMerLength) + 1; lCountBases++)
							{
								//  Get base at this position
								string strBase = strSequence.substr(lCountBases, 1);

								//  Compare the character, shift the sequence and append the appropriate base
								untSubSequence = untSubSequence << 2;
								if (strBase.compare("-") == 0)
								{
									ReportTimeStamp("[ProcessForeground16]", "ERROR:  Foreground Sequence Contains a Gap of Indeterminate Length");

									return false;
								}
								else if (strBase.compare("a") == 0)
									untSubSequence += m_untA;
								else if (strBase.compare("c") == 0)
									untSubSequence += m_untC;
								else if (strBase.compare("g") == 0)
									untSubSequence += m_untG;
								else if ((strBase.compare("t") == 0) || (strBase.compare("u") == 0))
									untSubSequence += m_untTU;
								//  Unknown character, this sequence is unsuitable to F-Distance analysis, mutate according to NCBI rules
								else if (bForegroundAllowUnknowns)
									untSubSequence += m_untA;
								else
								{
									ReportTimeStamp("[ProcessForeground16]", "ERROR:  Foreground Sequence Contains Unknown Characters and -background_allow_unknowns [-bau] is Not Set");

									return false;
								}

								//  Get mutation count
								nMutationCount = MutateSequence16(untSubSequence);
								if (nMutationCount >= 0)
								{
									//  Concatenate forward output
									if (nCountPass == 0)
										strForwardOutput += ConvertIntToString(nMutationCount);
									else
										strReverseOutput += ConvertIntToString(nMutationCount);
								}
								else
								{
									ReportTimeStamp("[ProcessForeground16]", "ERROR:  Foreground Sequence [" + ConvertSequenceToString16(untSubSequence, nNMerLength) + "] @ [" + ConvertLongToString(lCountBases) + "] Mutation Failed or Mutation Count Exceeds Limit [8]");

									return false;
								}
							}
						}

						//  Concatenate file text
						strOutputFileText = strForwardOutput;
						if (!strReverseOutput.empty())
							strOutputFileText += "\n" + strReverseOutput;

						//  Concatenate F-Distance table file text
						if (CompileFDistanceTableOutput(strAccession, strForwardOutput, strReverseOutput, strOutputTableEntry, false, false))
						{
							//  Write file text
							return WriteFileText(strOutputFilePathName, strOutputFileText);
						}
						else
						{
							ReportTimeStamp("[ProcessForeground16]", "ERROR:  F-Distance Score Compilation Failed");
						}
					}
					else
					{
						ReportTimeStamp("[ProcessForeground16]", "ERROR:  Foreground Sequence is Empty");
					}
				}
				else
				{
					ReportTimeStamp("[ProcessForeground16]", "ERROR:  Output File Path Name is Empty");
				}
			}
			else
			{
				ReportTimeStamp("[ProcessForeground16]", "ERROR:  Foreground Sequence Accession is Empty");
			}
		}
		else
		{
			ReportTimeStamp("[ProcessForeground16]", "ERROR:  Background Collection is Not Set");
		}
	}
	catch (exception ex)
	{
		cout << "ERROR [ProcessForeground16] Exception Code:  " << ex.what() << "\n";
	}

	return false;
}

////////////////////////////////////////////////////////////////////////////////
//
//  Markes a background sequence present in the background container
//         
////////////////////////////////////////////////////////////////////////////////
//
//  [uint16_t] untSequence:  sequenc to mark
//  [int] nNMerLength     :  nMer length to analyze
//                       :  returns true, if successful; else, false
//         
////////////////////////////////////////////////////////////////////////////////

bool MarkBackgroundSequence16(uint16_t& untSequence, int nNMerLength)
{
	try
	{
		//  If background array is set
		if (m_unaBackground16 != NULL)
		{
			//  Mark the first nMer in the background
			if ((untSequence >= 0) && (untSequence < UINT16_MAX))
			{
				omp_set_lock(&writelock16);
				m_unaBackground16[untSequence] = 1;
				omp_unset_lock(&writelock16);

				return true;
			}
			else if (untSequence == UINT32_MAX)
			{
				omp_set_lock(&writelock16);
				m_bBackgroundPolyTU16 = true;
				omp_unset_lock(&writelock16);

				return true;
			}
			else
			{
				ReportTimeStamp("[MarkBackgroundSequence16]", "ERROR:  Background Sequence [" + ConvertUnsignedInt16ToString(untSequence) + "] [" + ConvertSequenceToString16(untSequence, nNMerLength) + "] is Out of Range [0:" + ConvertUnsignedInt16ToString(UINT16_MAX) + "]");
			}
		}
		else
		{
			ReportTimeStamp("[MarkBackgroundSequence16]", "ERROR:  Background Container is Not Set");
		}
	}
	catch (exception ex)
	{
		cout << "ERROR [MarkBackgroundSequence16] Exception Code:  " << ex.what() << "\n";
	}

	return false;
}

////////////////////////////////////////////////////////////////////////////////
//
//  Mutates and marks a background sequence with unknown characters
//         
////////////////////////////////////////////////////////////////////////////////
//
//  [string] strSequence:  sequence to mutate
//  [int] nNMerLength   :  nMer length to analyze
//                     :  returns true, if successful; else, false
//         
////////////////////////////////////////////////////////////////////////////////

bool MutateAndMarkBackgroundSequence16(string& strSequence, int nNMerLength)
{
	//  Sequence numeric
	uint16_t untSequence = 0;
	//  Mutated sequence
	string strMutatedSequence = "";
	//  Base to mutate
	string strBase = "";
	//  Base max
	int nMaxBases = 0;
	//  Uknown type
	string strUnknownType = "";

	try
	{
		//  If background array is set
		if (m_unaBackground16 != NULL)
		{
			//  If the input sequence is not empty
			if (!strSequence.empty())
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

				//  If more unknowns are found, recurse
				if (strSequence.find_first_of('-') != string::npos)
				{
					ReportTimeStamp("[MutateAndMarkBackgroundSequence16]", "ERROR:  Background Sequence Contains a Gap of Indeterminate Length");

					return false;
				}
				else if ((strSequence.find_first_of('r') == string::npos) &&
					(strSequence.find_first_of('y') == string::npos) &&
					(strSequence.find_first_of('k') == string::npos) &&
					(strSequence.find_first_of('m') == string::npos) &&
					(strSequence.find_first_of('s') == string::npos) &&
					(strSequence.find_first_of('w') == string::npos) &&
					(strSequence.find_first_of('b') == string::npos) &&
					(strSequence.find_first_of('d') == string::npos) &&
					(strSequence.find_first_of('h') == string::npos) &&
					(strSequence.find_first_of('v') == string::npos) &&
					(strSequence.find_first_of('n') == string::npos))
				{
					untSequence = ConvertStringToSequence16(strSequence);
					return MarkBackgroundSequence16(untSequence, nNMerLength);
				}
				else
				{
					//  Set unknown type and max bases limit
					if (strSequence.find_first_of('r') != string::npos)
					{
						strUnknownType = "r";
						nMaxBases = 2;
					}
					else if (strSequence.find_first_of('y') != string::npos)
					{
						strUnknownType = "y";
						nMaxBases = 2;
					}
					else if (strSequence.find_first_of('k') != string::npos)
					{
						strUnknownType = "k";
						nMaxBases = 2;
					}
					else if (strSequence.find_first_of('m') != string::npos)
					{
						strUnknownType = "m";
						nMaxBases = 2;
					}
					else if (strSequence.find_first_of('s') != string::npos)
					{
						strUnknownType = "s";
						nMaxBases = 2;
					}
					else if (strSequence.find_first_of('w') != string::npos)
					{
						strUnknownType = "w";
						nMaxBases = 2;
					}
					else if (strSequence.find_first_of('b') != string::npos)
					{
						strUnknownType = "b";
						nMaxBases = 3;
					}
					else if (strSequence.find_first_of('d') != string::npos)
					{
						strUnknownType = "d";
						nMaxBases = 3;
					}
					else if (strSequence.find_first_of('h') != string::npos)
					{
						strUnknownType = "h";
						nMaxBases = 3;
					}
					else if (strSequence.find_first_of('v') != string::npos)
					{
						strUnknownType = "v";
						nMaxBases = 3;
					}
					else if (strSequence.find_first_of('n') != string::npos)
					{
						strUnknownType = "n";
						nMaxBases = 4;
					}

					//  Mutate one unknown character to all bases
					for (int nCountBases = 0; nCountBases < nMaxBases; nCountBases++)
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

						if (nCountBases == 0)
						{
							if ((strUnknownType == "r") || (strUnknownType == "w") || (strUnknownType == "d") || (strUnknownType == "h") || (strUnknownType == "m") || (strUnknownType == "v") || (strUnknownType == "n"))
								strBase = "a";
							else if ((strUnknownType == "y") || (strUnknownType == "s") || (strUnknownType == "b"))
								strBase = "c";
							else if (strUnknownType == "k")
								strBase = "g";
						}
						else if (nCountBases == 1)
						{
							if ((strUnknownType == "r") || (strUnknownType == "s") || (strUnknownType == "b") || (strUnknownType == "d"))
								strBase = "g";
							else if ((strUnknownType == "y") || (strUnknownType == "k") || (strUnknownType == "w"))
								strBase = "u";
							else if ((strUnknownType == "m") || (strUnknownType == "h") || (strUnknownType == "v") || (strUnknownType == "n"))
								strBase = "c";
						}
						else if (nCountBases == 2)
						{
							if ((strUnknownType == "b") || (strUnknownType == "d") || (strUnknownType == "h"))
								strBase = "u";
							else if ((strUnknownType == "v") || (strUnknownType == "n"))
								strBase = "g";
						}
						else if (nCountBases == 3)
							strBase = "u";

						strMutatedSequence = strSequence;
						strMutatedSequence.replace(strMutatedSequence.find_first_of(strUnknownType), 1, strBase);

						//  Else this is the last unknown character, mark the background array
						MutateAndMarkBackgroundSequence16(strMutatedSequence, nNMerLength);

						strBase = "";
					}
				}

				return true;
			}
			else
			{
				ReportTimeStamp("[MutateAndMarkBackgroundSequence16]", "ERROR:  Background Sequence is Empty");
			}
		}
		else
		{
			ReportTimeStamp("[MutateAndMarkBackgroundSequence16]", "ERROR:  Background Container is Not Set");
		}
	}
	catch (exception ex)
	{
		cout << "ERROR [MutateAndMarkBackgroundSequence16] Exception Code:  " << ex.what() << "\n";
	}

	return false;
}

////////////////////////////////////////////////////////////////////////////////
//
//  Processes the background genome; identifies background nMers present;
//         
////////////////////////////////////////////////////////////////////////////////
//
//  [CBase123_Catalog_Entry&] ceGet:  genome catalog entry
//  [string&] strSequence          :  sequence to process
//  [int] nNMerLength              :  nMer length to analyze
//  [bool] bBackgroundAllowUnknowns :  process background unknown chracters, if true
//                                :  returns true, if successful; else, false
//         
////////////////////////////////////////////////////////////////////////////////

bool ProcessBackground16(string& strSequence, int nNMerLength, bool bBackgroundAllowUnknowns)
{
	//  Sub-sequence
	string strSubSequence = "";
	//  The nMer sequence, binary
	uint16_t untSubSequence = 0;

	try
	{
		//  If background array is set
		if (m_unaBackground16 != NULL)
		{
			//  If the input sequence is not empty
			if (!strSequence.empty())
			{
				//  Reject sequences with gaps
				if (strSequence.find_first_of('-') != string::npos)
				{
					ReportTimeStamp("[ProcessBackground16]", "ERROR:  Background Sequence Contains a Gap of Indeterminate Length");

					return false;
				}
				//  Reject sequences with poly-n
				else if (strSequence.find("nnnnnnnn") != string::npos)
				{
					ReportTimeStamp("[ProcessBackground16]", "ERROR:  Input Sequence Contains Excessive Poly-n Bases");

					return false;
				}

				//  Get first nMer
				strSubSequence = strSequence.substr(0, nNMerLength);
				untSubSequence = ConvertStringToSequence16(strSubSequence);

				//  Process known bases
				if ((strSubSequence.find_first_of('r') == string::npos) &&
					(strSubSequence.find_first_of('y') == string::npos) &&
					(strSubSequence.find_first_of('k') == string::npos) &&
					(strSubSequence.find_first_of('m') == string::npos) &&
					(strSubSequence.find_first_of('s') == string::npos) &&
					(strSubSequence.find_first_of('w') == string::npos) &&
					(strSubSequence.find_first_of('b') == string::npos) &&
					(strSubSequence.find_first_of('d') == string::npos) &&
					(strSubSequence.find_first_of('h') == string::npos) &&
					(strSubSequence.find_first_of('v') == string::npos) &&
					(strSubSequence.find_first_of('n') == string::npos))
					MarkBackgroundSequence16(untSubSequence, nNMerLength);
				//  Process unknown bases
				else
				{
					if (bBackgroundAllowUnknowns)
					{
						if (!MutateAndMarkBackgroundSequence16(strSubSequence, nNMerLength))
						{
							ReportTimeStamp("[ProcessBackground16]", "ERROR:  Background Sequence Contains Unknown Characters and Unknown Characters Mutation Failed");

							return false;
						}
					}
					else
					{
						ReportTimeStamp("[ProcessBackground16]", "ERROR:  Background Sequence Contains Unknown Characters and -background_allow_unknowns [-bau] is Not Set");

						return false;
					}
				}

				//  Iterate subsequent nNMerLength characters to build remaing nMers
				for (long lCountBases = nNMerLength; lCountBases < (strSequence.length() - nNMerLength) + 1; lCountBases++)
				{
					//  Get base at this position
					string strBase = strSequence.substr(lCountBases, 1);

					//  Compare the character, shift the sequence and append the appropriate base
					untSubSequence = untSubSequence << 2;
					if (strBase.compare("-") == 0)
					{
						ReportTimeStamp("[ProcessBackground16]", "ERROR:  Background Sequence Contains a Gap of Indeterminate Length");

						return false;
					}
					else if (strBase.compare("a") == 0)
					{
						untSubSequence += m_untA;
						MarkBackgroundSequence16(untSubSequence, nNMerLength);
					}
					else if (strBase.compare("c") == 0)
					{
						untSubSequence += m_untC;
						MarkBackgroundSequence16(untSubSequence, nNMerLength);
					}
					else if (strBase.compare("g") == 0)
					{
						untSubSequence += m_untG;
						MarkBackgroundSequence16(untSubSequence, nNMerLength);
					}
					else if ((strBase.compare("t") == 0) || (strBase.compare("u") == 0))
					{
						untSubSequence += m_untTU;
						MarkBackgroundSequence16(untSubSequence, nNMerLength);
					}
					//  Unknown character, this sequence is unsuitable to F-Distance analysis, mutate according to NCBI rules
					else
					{
						if (bBackgroundAllowUnknowns)
						{
							strSubSequence = strSequence.substr(lCountBases - nNMerLength, nNMerLength);
							if (!MutateAndMarkBackgroundSequence16(strSubSequence, nNMerLength))
							{
								ReportTimeStamp("[ProcessBackground16]", "ERROR:  Background Sequence Contains Unknown Characters and Unknown Characters Mutation Failed");

								return false;
							}
						}
						else
						{
							ReportTimeStamp("[ProcessBackground16]", "ERROR:  Background Sequence Contains Unknown Characters and -background_allow_unknowns [-bau] is Not Set");

							return false;
						}
					}
				}

				return true;
			}
			else
			{
				ReportTimeStamp("[ProcessBackground16]", "ERROR:  Background Sequence is Empty");
			}
		}
		else
		{
			ReportTimeStamp("[ProcessBackground16]", "ERROR:  Background Container is Not Set");
		}
	}
	catch (exception ex)
	{
		cout << "ERROR [ProcessBackground16] Exception Code:  " << ex.what() << "\n";
	}

	return false;
}

////////////////////////////////////////////////////////////////////////////////
//
//  Processes the background/foreground genome file list to identify nMer(s) present;
//         
////////////////////////////////////////////////////////////////////////////////
//
//  [string] strInputListFilePathName     :  input file path name list
//  [string] strInputFilePathNameTransform:  input file path name transform (includes string replacements, see help)
//  [CBase123_Catalog&] b123Catalog       :  Base123 genome catalog to use
//  [bool] bBidirectional                 :  process bidirectional, if true
//  [int] nNMerLength                     :  nMer length to analyze
//  [bool] bBackground                    :  if true, process background; else, foreground
//  [bool] bAllowUnknowns                 :  process unknown chracters, if true
//  [string] strOutputFileNameSuffix      :  output file name suffix
//  [string] strErrorFilePathName         :  error file base name
//  [vector<string>&] vOutputTableEntries :  F-Distance table file text to concatenate
//  [int] nMaxProcs                       :  maximum processors for openMP
//                                       :  returns true, if successful; else, false
//         
////////////////////////////////////////////////////////////////////////////////

bool ProcessFDistanceList16(string strInputListFilePathName, string strInputFilePathNameTransform, CBase123_Catalog& b123Catalog, bool bBidirectional, int nNMerLength, bool bBackground, bool bAllowUnknowns,
	string strOutputFileNameSuffix, string strErrorFilePathName, vector<string>& vOutputTableEntries, int nMaxProcs)
{
	//  List lock
	omp_lock_t lockList;
	//  Input list file text
	string strInputListFileText = "";
	//  File path name vector<string>
	vector<string> vFilePathNames;
	//  Error file text
	vector<string> vErrorEntries;
	//  Error output file text
	string strErrorFileText = "";

	try
	{
		//  If input list file path name is not empty
		if (!strInputListFilePathName.empty())
		{
			//  If nMer length > 0
			if (nNMerLength > 0)
			{
				//  If background array is set
				if (m_unaBackground16 != NULL)
				{
					//  Get list file text
					if (GetFileText(strInputListFilePathName, strInputListFileText))
					{
						//  Split file path names
						SplitString(strInputListFileText, '\n', vFilePathNames);

						//  If vector contains file path names
						if (vFilePathNames.size() > 0)
						{
							//  Initialize output table vector if not backgrounbd
							if (!bBackground)
								vOutputTableEntries.resize(vFilePathNames.size());

							//  Initialize error file vector
							vErrorEntries.resize(vFilePathNames.size());

							//  Initialize time stamp lock
							omp_init_lock(&lockList);

							//  Declare team size
							#pragma omp parallel shared(m_unaBackground16, vOutputTableEntries) num_threads(nMaxProcs)
							{
								#pragma omp for
								for (long lCount = 0; lCount < vFilePathNames.size(); lCount++)
								{
									//  Test max procs
									if (lCount == 0)
									{
										omp_set_lock(&lockList);
										ReportTimeStamp("[ProcessFDistanceList16]", "NOTE:  Thread Count = " + ConvertIntToString(omp_get_num_threads()));
										omp_unset_lock(&lockList);
									}

									//  Update for timestamp every 10,000 files
									if (lCount % 10000 == 0)
									{
										omp_set_lock(&lockList);
										ReportTimeStamp("[ProcessFDistanceList16]", "NOTE:  Processing Entry [" + ConvertLongToString(lCount) + "] [" + vFilePathNames[lCount] + "]");
										omp_unset_lock(&lockList);
									}


									//  If the file  name is not empty
									if (!vFilePathNames[lCount].empty())
									{
										//  Sequence file text
										string strSequenceFileText = "";
										//  Working file path name
										string strWorkingFilePathName = "";

										//  If input file path name transform is not empty
										if (!strInputFilePathNameTransform.empty())
											strWorkingFilePathName = TransformFilePathName(vFilePathNames[lCount], strInputFilePathNameTransform, "");
										else
											strWorkingFilePathName = vFilePathNames[lCount];

										//  Get sequence file text
										if (GetFileText(strWorkingFilePathName, strSequenceFileText))
										{
											//  Accession
											string strAccession = "";

											//  Get Accession
											strAccession = GetAccessionFromFileHeader(strSequenceFileText);

											if (!strAccession.empty())
											{
												//  Catalog entry
												CBase123_Catalog_Entry ceGet;

												//vbContinue = b123Catalog.GetEntryByAccession(strAccession, ceGet);
												if (b123Catalog.GetEntryByAccession(strAccession, ceGet))
												{
													//  Sequence
													string strForward = "";

													//  Get sequence
													strForward = GetSequenceFromFAFile(strSequenceFileText);

													//  Process forward sequence
													if (!strForward.empty())
													{
														// If sequence is circular
														if (ceGet.GetStrandednessType() == "c")
														{
															//  Circularize
															strForward += strForward.substr(0, nNMerLength - 1);
														}

														//  Process background sequence
														if (bBackground)
														{
															//  Process background forward
															//vbContinue = ProcessBackground16(strForward, nNMerLength, bAllowUnknowns);
															if (ProcessBackground16(strForward, nNMerLength, bAllowUnknowns))
															{
																//  If bidirectional processing required
																if (bBidirectional)
																{
																	//  Reverse compliment
																	string strReverse = "";

																	//  Get reverse compliment
																	strReverse = ConvertToReverseCompliment(strForward);

																	//  Process background reverse compliment
																	//vbContinue = ProcessBackground16(strReverse, nNMerLength, bAllowUnknowns);
																	if (!ProcessBackground16(strReverse, nNMerLength, bAllowUnknowns))
																	{
																		vErrorEntries[lCount] = strWorkingFilePathName + "~Background (Reverse) Analysis Failed\n";

																		omp_set_lock(&lockList);
																		ReportTimeStamp("[ProcessFDistanceList16]", "ERROR:  Background Reverse Sequence [" + strWorkingFilePathName + "] Process Failed");
																		omp_unset_lock(&lockList);
																	}
																}
															}
															else
															{
																vErrorEntries[lCount] = strWorkingFilePathName + "~Background (Forward) Analysis Failed\n";

																omp_set_lock(&lockList);
																ReportTimeStamp("[ProcessFDistanceList16]", "ERROR:  Background Forward Sequence [" + strWorkingFilePathName + "] Process Failed");
																omp_unset_lock(&lockList);
															}
														}
														//  Process foreground sequence
														else
														{
															//  Output file path name
															string strOutputFilePathName = "";
															//  Path delimiter
															char chrPathDelimiter = strPathDelimiter;

															//  Get base path name and concatenate output file path name
															if(!strOutputFileNameSuffix.empty())
																strOutputFilePathName = GetBasePath(strWorkingFilePathName) + chrPathDelimiter + GetFileNameExceptLastExtension(strWorkingFilePathName) + "." + strOutputFileNameSuffix + ".fdist";
															else
																strOutputFilePathName = GetBasePath(strWorkingFilePathName) + chrPathDelimiter + GetFileNameExceptLastExtension(strWorkingFilePathName) + ".fdist";

															if (!ProcessForeground16(strOutputFilePathName, strAccession, strForward, bBidirectional, bAllowUnknowns, nNMerLength, vOutputTableEntries[lCount]))
															{
																vErrorEntries[lCount] = strWorkingFilePathName + "~Foreground Analysis Failed\n";

																omp_set_lock(&lockList);
																ReportTimeStamp("[ProcessFDistanceList16]", "ERROR:  Foreground Sequence [" + strWorkingFilePathName + "] Process Failed");
																omp_unset_lock(&lockList);
															}
														}
													}
													else
													{
														vErrorEntries[lCount] = strWorkingFilePathName + "~Empty Sequence\n";

														omp_set_lock(&lockList);
														ReportTimeStamp("[ProcessFDistanceList16]", "ERROR:  Input Sequence File [" + strWorkingFilePathName + "] Sequence is Empty");
														omp_unset_lock(&lockList);
													}
												}
												else
												{
													vErrorEntries[lCount] = strWorkingFilePathName + "~Catalog Accession Search Failed\n";

													omp_set_lock(&lockList);
													ReportTimeStamp("[ProcessFDistanceList16]", "ERROR:  Input Sequence File [" + strWorkingFilePathName + "] Catalog Accession [" + strAccession + "] Search Failed");
													omp_unset_lock(&lockList);
												}

												//  Clear accession
												strAccession = "";
											}
											else
											{
												vErrorEntries[lCount] = strWorkingFilePathName + "~Empty Accession\n";

												omp_set_lock(&lockList);
												ReportTimeStamp("[ProcessFDistanceList16]", "ERROR:  Input Sequence File [" + strWorkingFilePathName + "] Accession is Empty");
												omp_unset_lock(&lockList);
											}
										}
										else
										{
											vErrorEntries[lCount] = strWorkingFilePathName + "~File Open Failed\n";

											omp_set_lock(&lockList);
											ReportTimeStamp("[ProcessFDistanceList16]", "ERROR:  Input Sequence File [" + strWorkingFilePathName + "] Open Failed");
											omp_unset_lock(&lockList);
										}
									}
									//  Report no error, empty lines should not be present
								}
							}

							//  Destroy time stamp lock
							omp_destroy_lock(&lockList);

							//  Write error file
							if (!strErrorFilePathName.empty())
							{
								//  Add header
								strErrorFileText = "File Path Name~Error\n";

								//  Iterate error entries and concatenate error file text
								for (long lCount = 0; lCount < vErrorEntries.size(); lCount++)
								{
									//  If file error entry is not empty, concatenate error file text
									if (!vErrorEntries[lCount].empty())
										strErrorFileText += vErrorEntries[lCount];
								}

								//  Write error file
								WriteFileText(strErrorFilePathName, strErrorFileText);
							}

							vFilePathNames.clear();
							vErrorEntries.clear();

							return true;
						}
						else
						{
							ReportTimeStamp("[ProcessFDistanceList16]", "ERROR:  Input File Path Name List [" + strInputListFilePathName + "] Text is Empty");
						}
					}
					else
					{
						ReportTimeStamp("[ProcessFDistanceList16]", "ERROR:  Input File Path Name List [" + strInputListFilePathName + "] Open Failed");
					}
				}
				else
				{
					ReportTimeStamp("[ProcessFDistanceList16]", "ERROR:  Background Container is Not Set");
				}
			}
			else
			{
				ReportTimeStamp("[ProcessFDistanceList16]", "ERROR:  nMer Length is Not Properly Set:  Should be 8 or 16");
			}
		}
		else
		{
			ReportTimeStamp("[ProcessFDistanceList16]", "ERROR:  Input File Path Name List is Empty");
		}
	}
	catch (exception ex)
	{
		cout << "ERROR [ProcessFDistanceList16] Exception Code:  " << ex.what() << "\n";
	}

	return false;
}

////////////////////////////////////////////////////////////////////////////////
//
//  Initializes the background array
//         
////////////////////////////////////////////////////////////////////////////////
//
//  :  returns true, if successful; else, false
//         
////////////////////////////////////////////////////////////////////////////////

bool InitializeBackground16()
{
	try
	{
		//  Initialize the background array
		if (m_unaBackground16 == NULL)
			m_unaBackground16 = new uint16_t[UINT16_MAX];

		//  If background array is set
		if (m_unaBackground16 != NULL)
		{
			//  Initialize background to 0
			for (uint16_t untCount = 0; untCount < UINT16_MAX; untCount++)
				m_unaBackground16[untCount] = 0;

			return true;
		}
		else
		{
			ReportTimeStamp("[InitializeBackground16]", "ERROR:  16-mer Background Initialization Failed");
		}
	}
	catch (exception ex)
	{
		cout << "ERROR [InitializeBackground16] Exception Code:  " << ex.what() << "\n";
	}

	return false;
}

////////////////////////////////////////////////////////////////////////////////
//
//  Destroys the background array
//         
////////////////////////////////////////////////////////////////////////////////
//
//  :  returns true, if successful; else, false
//         
////////////////////////////////////////////////////////////////////////////////

bool DestroyBackground16()
{
	try
	{
		//  If background array is set
		if (m_unaBackground16 != NULL)
			delete[] m_unaBackground16;

		return true;
	}
	catch (exception ex)
	{
		cout << "ERROR [DestroyBackground16] Exception Code:  " << ex.what() << "\n";
	}

	return false;
}

////////////////////////////////////////////////////////////////////////////////
//
//  Initializes the write lock
//         
////////////////////////////////////////////////////////////////////////////////
//
//  :  returns true, if successful; else, false
//         
////////////////////////////////////////////////////////////////////////////////

bool InitializeWriteLock16()
{
	try
	{
		//  Initialize write lock
		omp_init_lock(&writelock16);

		return true;
	}
	catch (exception ex)
	{
		cout << "ERROR [InitializeWriteLock16] Exception Code:  " << ex.what() << "\n";
	}

	return false;
}

////////////////////////////////////////////////////////////////////////////////
//
//  Destroys the write lock
//         
////////////////////////////////////////////////////////////////////////////////
//
//  :  returns true, if successful; else, false
//         
////////////////////////////////////////////////////////////////////////////////

bool DestroyWriteLock16()
{
	try
	{
		//  Destroy write lock
		omp_destroy_lock(&writelock16);

		return true;
	}
	catch (exception ex)
	{
		cout << "ERROR [DestroyWriteLock16] Exception Code:  " << ex.what() << "\n";
	}

	return false;
}

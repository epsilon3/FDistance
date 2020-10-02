// F_Dist_R.h : Declare member variables

////////////////////////////////////////////////////////////////////////////////
//
//  F_Dist_R wrapper main application (header):
//
//  Developed by Stephen Donald Huff, PhD (Stephen.Huff.3@us.af.mil)
//       and William Chan (williamchan@roadrunner.com)
//  Biological Informatics Group, RHDJ, 711HPW, United States Air Force Research Laboratory
//  18 July 2017
//  (All Rights Reserved)
//
////////////////////////////////////////////////////////////////////////////////

#pragma once

using namespace std;

#include <string>
#include <vector>
#include <iostream>

//  Base123 members

#ifdef _WIN64
#define strPathDelimiter '\\';
#else
#define strPathDelimiter '/';
#endif

//  CDS structure
struct structCDS
{
	string strIsComplement;
	string strNameID;
	long lStart;
	long lStop;
	string strCompleteness;
};

//  Adenine nucleotide value
const unsigned int m_untA = 0b00;
//  Thymine/uracil nucleotide value
const unsigned int m_untC = 0b01;
//  Guanine nucleotide value
const unsigned int m_untG = 0b10;
//  Cytosine nucleotide value
const unsigned int m_untTU = 0b11;
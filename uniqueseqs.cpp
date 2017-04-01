#include <iostream>
#include <string>
#include <vector>
#include <fstream>
#include <algorithm>

using namespace std ;

// Comparison; not case sensitive.
bool compareNoCase (string first, string second)
{
  int i=0;
  while ((i < first.length()) && (i < second.length()))
  {
    if (tolower (first[i]) < tolower (second[i])) return true;
    else if (tolower (first[i]) > tolower (second[i])) return false;
    i++;
  }

  if (first.length() < second.length()) return true;
  else return false;
}


int main () 
{

    ///////////////////////////////////////////////////////////////////////////////////
    // Declaration of variables

    ifstream ifile ;

    ofstream ofile ;

	string line ;

	vector <string> dnaseq ;

    ///////////////////////////////////////////////////////////////////////////////////
    // Opens all necessary files.
	
    ifile.open("out");

    ofile.open("unique.txt") ;

    ///////////////////////////////////////////////////////////////////////////////////
    // Takes all data from in file, reads it line by line and puts each line into a
    // vector element.
    int size = 209 ;
    // Modify this size according to what type of data you want to extract from the file.
    while( ifile >> line )
    {
        if( line > size)
        {
            //cout << line << endl ;
            dnaseq.push_back(line) ;
        }
    }

    ///////////////////////////////////////////////////////////////////////////////////
    // Sorts all vector elements alphabetically, puts duplicate elements at the end of
    //the vector list and erases them.

    sort( dnaseq.begin(), dnaseq.end(), compareNoCase );

    dnaseq.erase( unique( dnaseq.begin(), dnaseq.end() ) , dnaseq.end());

    ///////////////////////////////////////////////////////////////////////////////////
    // Outputs vector contents into a new file.

    for( int i = 0 ; i < dnaseq.size() ; i++ )
    {
        ofile << dnaseq[i] << endl ;
    }

    ///////////////////////////////////////////////////////////////////////////////////
    // Closes all necessary files.

    ifile.close() ;

    ofile.close() ;

    ///////////////////////////////////////////////////////////////////////////////////

	return 0 ;

}

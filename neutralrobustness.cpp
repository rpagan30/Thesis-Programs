#include<iostream>
#include<cmath>
#include<stdio.h>
#include<stdlib.h>
#include<fstream>
#include<string>
#include<iomanip>
#include<vector>
//old random number generator:
#include "randomc.h" 
#include "mersenne.cpp"
#include "userintf.cpp"

//generates a Poisson distribution:
#include "stocc.h"
#include "stoc1.cpp"
#include<time.h>

//new random number  generator
#include "./Mersenne-1.1/MersenneTwister.h"

//the previous code used a random number  generator (mersenne.cpp) that was very slow and would not generate the appropriate numbers initially for the nucleotide positions that were mutated. This was also true on 32 bit machines. This was fixed by using MersenneTwister.h 

///program structure : a pdb file has its energy refined and calculated intially. Then this amino acid sequence is used to populate a population. Mutations are introduced and REFINED energies calculated scrwl3 is the key to this program as it will replace the heavy atoms of a mutated residue with new coordinates for the new residue


// to run the program : ./neutralrobustness

using namespace std;

#define GEN 500
#define POPN 1000
#define INIT_RANDOM 0
#define RANDOM_GENERATOR TRandomMersenne

MTRand mtrand1;

char AA[20] = {'A','E','Q','D','N','L','G','K','S','V','R','T','P','I','M','F','Y','C','W','H'};


int seed = time(0);
TRandomMersenne rg(seed);
ofstream robustout("robustout");
ofstream energyout("energyout");
ofstream outfile("out");
/////////////////the following generates the energy for the protein///////////////////////////

float energy(vector<float>& x, vector<float>& y, vector<float>& z,vector<int>& num, int n, vector<char>& seq, float contactenergy[20][20])
{
int aacontacts[1000][1000] = {0};
 float dist;
 float eng = 0.0;
 int f=0,m=0;

	for(int i=0;i<n;i++)
	   { 
		for(int j=0;j<n;j++)
		   {
		   dist = (x[i]-x[j])*(x[i]-x[j])+(y[i]-y[j])*(y[i]-y[j])+(z[i]-z[j])*(z[i]-z[j]);
		    if(dist<=20.25 && abs(num[i]-num[j])>=2) aacontacts[num[i]][num[j]] = 1;
		    //else aacontacts[num[i]][num[j]] = 0;
         
		   }
	    }	

/////extract the amino acid sequence from num and seq////////////

char aminoacidseq[70];
int p=0;
     
for(int q=0;q<num.size();q++)
   {
    if(num[q]!=num[q+1]|| num.size() < q + 2 )
      {aminoacidseq[p] = seq[q];
       p++;}
   }
//////////calculate energy//////////////////////////////////////////////////////////
	for(int i=1; i<71;i++)
	{
		for(int j=i;j<71;j++) /*ale was here*/
		{
			if(aacontacts[i][j] ==1)

			{
			//cout<<i<<" "<<j<<endl;
				for(int l=0;l<20;l++)
				{
					if(AA[l] == aminoacidseq[i-1])
						m = l;
					if(AA[l] == aminoacidseq[j-1])
						f = l;					
				}
				eng=eng+contactenergy[m][f];
				//cout<<contactenergy[m][f]<<endl;
				//cout<<AA[m]<<" "<<AA[f]<<endl;
				m=0;
				f=0;
				contactenergy[m][f] = 0.0;
				//cout<<eng<<endl;
			}
		}
	}
        //eng = (eng/2); /* ale was here */
	//cout<<eng<<"\n"<<endl;
	return eng;
}


/////////////the following translates the DNA sequence/////////////////////////////////////////////////////////////////////////

char *translate(string codon)
{

	if      (codon == "TCA")    { return "s";}   
    	else if (codon == "TCC")    { return "s";}   
	else if (codon == "TCG")    { return "s";}  
    	else if (codon == "TCT")    { return "s";}    
    	else if (codon == "TTC")    { return "f";}
    	else if (codon == "TTT")    { return "f";}    
    	else if (codon == "TTA")    { return "l";} 
    	else if (codon == "TTG")    { return "l";}   
    	else if (codon == "TAC")    { return "y";} 
    	else if (codon == "TAT")    { return "y";}  
    	else if (codon == "TAA")    { return "_";} 
    	else if (codon == "TAG")    { return "_";} 
    	else if (codon == "TGC")    { return "c";}  
    	else if (codon == "TGT")    { return "c";}  
    	else if (codon == "TGA")    { return "_";}  
    	else if (codon == "TGG")    { return "w";} 
    	else if (codon == "CTA")    { return "l";}  
    	else if (codon == "CTC")    { return "l";}  
    	else if (codon == "CTG")    { return "l";}  
    	else if (codon == "CTT")    { return "l";}  
    	else if (codon == "CCA")    { return "p";}   
    	else if (codon == "CCC")    { return "p";}   
    	else if (codon == "CCG")    { return "p";}   
    	else if (codon == "CCT")    { return "p";}  
    	else if (codon == "CAC")    { return "h";}  
    	else if (codon == "CAT")    { return "h";}  
    	else if (codon == "CAA")    { return "q";}   
    	else if (codon == "CAG")    { return "q";}   
    	else if (codon == "CGA")    { return "r";}   
    	else if (codon == "CGC")    { return "r";}  
    	else if (codon == "CGG")    { return "r";}   
    	else if (codon == "CGT")    { return "r";}   
    	else if (codon == "ATA")    { return "i";}  
    	else if (codon == "ATC")    { return "i";}   
    	else if (codon == "ATT")    { return "i";}   
    	else if (codon == "ATG")    { return "m";}   
    	else if (codon == "ACA")    { return "t";}    
    	else if (codon == "ACC")    { return "t";}    
    	else if (codon == "ACG")    { return "t";}   
    	else if (codon == "ACT")    { return "t";}   
    	else if (codon == "AAC")    { return "n";}    
    	else if (codon == "AAT")    { return "n";}    
    	else if (codon == "AAA")    { return "k";}   
    	else if (codon == "AAG")    { return "k";}   
    	else if (codon == "AGC")    { return "s";}   
    	else if (codon == "AGT")    { return "s";}   
    	else if (codon == "AGA")    { return "r";}   
    	else if (codon == "AGG")    { return "r";}  
    	else if (codon == "GTA")    { return "v";}   
    	else if (codon == "GTC")    { return "v";}   
    	else if (codon == "GTG")    { return "v";}    
    	else if (codon == "GTT")    { return "v";}   
    	else if (codon == "GCA")    { return "a";}    
    	else if (codon == "GCC")    { return "a";}   
    	else if (codon == "GCG")    { return "a";}    
    	else if (codon == "GCT")    { return "a";}    
    	else if (codon == "GAC")    { return "d";} 
    	else if (codon == "GAT")    { return "d";}  
    	else if (codon == "GAA")    { return "e";}  
    	else if (codon == "GAG")    { return "e";}    
    	else if (codon == "GGA")    { return "g";}   
    	else if (codon == "GGC")    { return "g";}    
    	else if (codon == "GGG")    { return "g";}    
    	else if (codon == "GGT")    { return "g";}   
    	else
	{
	        cout<<"Unrecognized codon: "<<codon<<endl;
    	}

}

/////////////////////////////the following is the mutation procedure////////////////////////////////////////////////////////////


string Mutate(string d, int m)
{
int nuc, nu, trans;
    if(m!=0)
      {for(int i=0; i<m; i++)
        {
           //nuc = rg.IRandom(0,d.length()-1); //cout<<"random nuc = "<<nuc<<endl;
           nuc = mtrand1.randInt(d.length()-1);
           //nu = rg.IRandom(0,2); //cout<<"new nuc = "<<nu<<endl;
           nu = mtrand1.randInt(2);
           //trans = rg.IRandom(0,1); //cout<<"trans = "<<trans<<endl;
           trans = mtrand1.randInt(2);

           if(d[nuc] == 'T')
             {
              if(nu == 0 || nu == 1)
                {d[nuc] = 'C';}
              else {if(trans == 0)
                       {d[nuc] = 'A';}
                    else d[nuc] = 'G';
                   }              
             }
           else if(d[nuc] == 'A')
             {
              if(nu == 0 || nu == 1)
                {d[nuc] = 'G';}
              else {if(trans == 0)
                       {d[nuc] = 'T';}
                    else d[nuc] = 'C';
                   }
             }
           else if(d[nuc] == 'C')
             {
             if(nu == 0 || nu == 1)
                {d[nuc] = 'T';}
              else {if(trans == 0)
                       {d[nuc] = 'A';}
                    else d[nuc] = 'G';
                   }
             }
           else if(d[nuc] == 'G')
             {
             if(nu == 0 || nu == 1)
                {d[nuc] = 'A';}
              else {if(trans == 0)
                       {d[nuc] = 'T';}
                    else d[nuc] = 'C';
                   }
             }
        }
      }

return d;
}

//////////////////////the following calculates the robustness of the most common sequence in the population /////////////////////////

float robustness(float contactenergy[20][20])
{
std::vector<int> num;
std::vector<float> X;
std::vector<float> Y;
std::vector<float> Z;
std::vector<char> aaseq;
string line;

string aainput, dnainput, dnainput2, cod, dnanew, mutatedcodon;
system("./mostcommonDNA.pl popndnaseqs > mostcommonout");

ifstream comdnafile;
comdnafile.open("mostcommonout");
comdnafile>>dnainput;
//robustout<<dnainput<<endl;

////////////////////calculate the energy of the most common DNA sequence//////////////////////////////////////////////////////

string protein, dnaseq, codon;
ofstream aafile("aa2");
//system("rm perlinput");
ofstream perlinput("perlinput");

dnaseq = "";protein = "";
dnaseq = dnainput;
    for(int z=0; z<(dnaseq.length()-2);z += 3)
       {
	codon = dnaseq.substr(z,3);
       	if(codon == "TAA"|| codon == "TAG" ||codon =="TGA")
	  {
	  break;
	  }
        else protein += translate(codon); 
   	}

aafile<<protein<<endl;


//system("./scwrl3_lin/scwrl3 -i 1EGL.pdb -o refined.pdb -s aa2 > logfile");
system("./scwrl4_lin/Scwrl4 -i 1EGL.pdb -o refined.pdb -s aa2 -h -v > logfile");
float origenergy;
ifstream pdbfile;
pdbfile.open("refined.pdb");

  if (pdbfile.is_open())
  {
    while (! pdbfile.eof())
    {
      getline (pdbfile,line);
      int Epos = line.find("ATOM",0);
      int Epos2 = line.find("A",21);

      
       if(Epos2 != string::npos && Epos != string::npos)
       {
      // cout<<line<<endl;
       string str = line.substr(17,3);
      
       if(str == "ALA"){/*cout<<"A"<<endl;*/ aaseq.push_back('A');}
       if(str == "ARG"){/*cout<<"R"<<endl;*/ aaseq.push_back('R');}
       if(str == "ASN"){/*cout<<"N"<<endl;*/ aaseq.push_back('N');}
       if(str == "ASP"){/*cout<<"D"<<endl;*/ aaseq.push_back('D');}
       if(str == "CYS"){/*cout<<"C"<<endl;*/ aaseq.push_back('C');}
       if(str == "GLN"){/*cout<<"Q"<<endl;*/ aaseq.push_back('Q');}
       if(str == "GLU"){/*cout<<"E"<<endl;*/ aaseq.push_back('E');}
       if(str == "GLY"){/*cout<<"G"<<endl;*/ aaseq.push_back('G');}
       if(str == "HIS"){/*cout<<"H"<<endl;*/ aaseq.push_back('H');}
       if(str == "ILE"){/*cout<<"I"<<endl;*/ aaseq.push_back('I');}
       if(str == "LEU"){/*cout<<"L"<<endl;*/ aaseq.push_back('L');}
       if(str == "LYS"){/*cout<<"K"<<endl;*/ aaseq.push_back('K');}
       if(str == "MET"){/*cout<<"M"<<endl;*/ aaseq.push_back('M');}
       if(str == "PHE"){/*cout<<"F"<<endl;*/ aaseq.push_back('F');}
       if(str == "PRO"){/*cout<<"P"<<endl;*/ aaseq.push_back('P');}
       if(str == "SER"){/*cout<<"S"<<endl;*/ aaseq.push_back('S');}
       if(str == "THR"){/*cout<<"T"<<endl;*/ aaseq.push_back('T');}
       if(str == "TRP"){/*cout<<"W"<<endl;*/ aaseq.push_back('W');}
       if(str == "TYR"){/*cout<<"Y"<<endl;*/ aaseq.push_back('Y');}
       if(str == "VAL"){/*cout<<"V"<<endl;*/ aaseq.push_back('V');}

 	       string seqnum = line.substr(23, 3);
   	       for (int i = 0; i < seqnum.length(); i++)
     		     {
    	          if (seqnum[i] == ' ')
        	      {
         	     seqnum.replace(i,1,"");
         	     i = i - 1;
       	    		   } 
       		     }
      		int seqnumint = atoi(seqnum.c_str());
     		num.push_back(seqnumint);

      	 	string x = line.substr(31, 7);
         	 for (int i = 0; i < x.length(); i++)
        	  {
          	    if (x[i] == ' ')
          	    {
          	    x.replace(i,1,"");
          	    i = i - 1;
          	    } 
       	          }
      		float xint = atof(x.c_str());
		X.push_back(xint);
 
  	       string y = line.substr(39,7);
   	       for (int i = 0; i < y.length(); i++)
   	       {
      	        if (y[i] == ' ')
      	        {
       	         y.replace(i,1,"");
        	 i = i - 1;
         	 } 
               }
    	       float yint = atof(y.c_str());		Y.push_back(yint);


      		string z = line.substr(47,7);
         	 for (int i = 0; i < z.length(); i++)
        	  {
        	     if (z[i] == ' ')
          	     {
         	     z.replace(i,1,"");
          	     i = i - 1;
          	     } 
                  }
                float zint = atof(z.c_str());
 		Z.push_back(zint);
          }
      }
    pdbfile.close();
    remove("refined.pdb");
    remove("aa2");
   }


origenergy = energy(X,Y,Z,num,aaseq.size(),aaseq, contactenergy);
energyout<<origenergy<<endl;
perlinput<<"initial refined energy = "<<origenergy<<endl;
perlinput<<"DNA Seq : "<<dnainput<<endl;
perlinput<<"AA Seq : "<<protein<<endl;

//////////put the point mutated dna sequences into an array////////////////////////////////////////

std::vector<string> dnalist;
mutatedcodon ="";

for(int i=0;i<(dnainput.length()-2);i += 3)
   {
    cod = dnainput.substr(i,3);
    string one = cod.substr(0,1);
    string two = cod.substr(1,1); 
    string three = cod.substr(2,1);

//////////first codon position:

	if(cod[0] == 'A')
 	  {
          mutatedcodon = "T" + two + three; 
          dnainput2 = dnainput;
          dnanew = dnainput2.replace(i,3,mutatedcodon);
  	  dnalist.push_back(dnanew);
          dnanew = "";
          mutatedcodon = "";
          
          mutatedcodon = "C" + two + three; 
          dnainput2 = dnainput;
          dnanew = dnainput2.replace(i,3,mutatedcodon);
  	  dnalist.push_back(dnanew);
          dnanew = "";
          mutatedcodon = "";
  	  
          mutatedcodon = "G" + two + three; 
          dnainput2 = dnainput;
          dnanew = dnainput2.replace(i,3,mutatedcodon);
  	  dnalist.push_back(dnanew);
          dnanew = "";
          mutatedcodon = "";
    	  }

	if(cod[0] == 'T')
 	  {
          mutatedcodon = "A" + two + three; 
          dnainput2 = dnainput;
          dnanew = dnainput2.replace(i,3,mutatedcodon);
  	  dnalist.push_back(dnanew);
          dnanew = "";
          mutatedcodon = "";
          
          mutatedcodon = "C" + two + three; 
          dnainput2 = dnainput;
          dnanew = dnainput2.replace(i,3,mutatedcodon);
  	  dnalist.push_back(dnanew);
          dnanew = "";          mutatedcodon = "";
  	  
          mutatedcodon = "G" + two + three; 
          dnainput2 = dnainput;
          dnanew = dnainput2.replace(i,3,mutatedcodon);
  	  dnalist.push_back(dnanew);
          dnanew = "";
          mutatedcodon = "";
    	  }

    	if(cod[0] == 'G')
 	  {          mutatedcodon = "A" + two + three; 
          dnainput2 = dnainput;
          dnanew = dnainput2.replace(i,3,mutatedcodon);
  	  dnalist.push_back(dnanew);
          dnanew = "";
          mutatedcodon = "";
          
          mutatedcodon = "C" + two + three; 
          dnainput2 = dnainput;
          dnanew = dnainput2.replace(i,3,mutatedcodon);
  	  dnalist.push_back(dnanew);
          dnanew = "";
          mutatedcodon = "";
  	  
          mutatedcodon = "T" + two + three;           dnainput2 = dnainput;
          dnanew = dnainput2.replace(i,3,mutatedcodon);
  	  dnalist.push_back(dnanew);
          dnanew = "";
          mutatedcodon = "";
    	  }

   	if(cod[0] == 'C')
 	  {
          mutatedcodon = "A" + two + three; 
          dnainput2 = dnainput;
          dnanew = dnainput2.replace(i,3,mutatedcodon);
  	  dnalist.push_back(dnanew);
          dnanew = "";
          mutatedcodon = "";
          
          mutatedcodon = "G" + two + three; 
          dnainput2 = dnainput;
          dnanew = dnainput2.replace(i,3,mutatedcodon);
  	  dnalist.push_back(dnanew);
          dnanew = "";
          mutatedcodon = "";
  	  
          mutatedcodon = "T" + two + three; 
          dnainput2 = dnainput;
          dnanew = dnainput2.replace(i,3,mutatedcodon);
  	  dnalist.push_back(dnanew);
          dnanew = "";
          mutatedcodon = "";
    	  }


///second codon position:

	if(cod[1] == 'A')
 	  {
          mutatedcodon = one + "T" + three; 
          dnainput2 = dnainput;
          dnanew = dnainput2.replace(i,3,mutatedcodon);
  	  dnalist.push_back(dnanew);          dnanew = "";
          mutatedcodon = "";
          
          mutatedcodon = one +"C" + three; 
          dnainput2 = dnainput;
          dnanew = dnainput2.replace(i,3,mutatedcodon);
  	  dnalist.push_back(dnanew);
          dnanew = "";
          mutatedcodon = "";
  	  
          mutatedcodon = one + "G" + three; 
          dnainput2 = dnainput;
          dnanew = dnainput2.replace(i,3,mutatedcodon);
  	  dnalist.push_back(dnanew);
          dnanew = "";
          mutatedcodon = "";
    	  }

	if(cod[1] == 'T')
 	  {
          mutatedcodon = one + "A" + three; 
          dnainput2 = dnainput;
          dnanew = dnainput2.replace(i,3,mutatedcodon);
  	  dnalist.push_back(dnanew);
          dnanew = "";
          mutatedcodon = "";
          
          mutatedcodon = one + "C" + three; 
          dnainput2 = dnainput;
          dnanew = dnainput2.replace(i,3,mutatedcodon);
  	  dnalist.push_back(dnanew);
          dnanew = "";
          mutatedcodon = "";
  	  
          mutatedcodon = one + "G" + three; 
          dnainput2 = dnainput;
          dnanew = dnainput2.replace(i,3,mutatedcodon);
  	  dnalist.push_back(dnanew);
          dnanew = "";
          mutatedcodon = "";
    	  }

    	if(cod[1] == 'G')
 	  {
          mutatedcodon = one + "A" + three; 
          dnainput2 = dnainput;
          dnanew = dnainput2.replace(i,3,mutatedcodon);
  	  dnalist.push_back(dnanew);
          dnanew = "";
          mutatedcodon = "";
          
          mutatedcodon = one + "C" + three; 
          dnainput2 = dnainput;
          dnanew = dnainput2.replace(i,3,mutatedcodon);
  	  dnalist.push_back(dnanew);
          dnanew = "";
          mutatedcodon = "";
  	  
          mutatedcodon = one + "T" + three; 
          dnainput2 = dnainput;
          dnanew = dnainput2.replace(i,3,mutatedcodon);
  	  dnalist.push_back(dnanew);
          dnanew = "";
          mutatedcodon = "";
    	  }

   	if(cod[1] == 'C')
 	  {
          mutatedcodon = one + "A" + three; 
          dnainput2 = dnainput;
          dnanew = dnainput2.replace(i,3,mutatedcodon);
  	  dnalist.push_back(dnanew);
          dnanew = "";
          mutatedcodon = "";
          
          mutatedcodon = one + "G" + three;           dnainput2 = dnainput;
          dnanew = dnainput2.replace(i,3,mutatedcodon);
  	  dnalist.push_back(dnanew);
          dnanew = "";
          mutatedcodon = "";
  	  
          mutatedcodon = one + "T" + three; 
          dnainput2 = dnainput;
          dnanew = dnainput2.replace(i,3,mutatedcodon);
  	  dnalist.push_back(dnanew);
          dnanew = "";
          mutatedcodon = "";
    	  }

/////third codon position

	if(cod[2] == 'A')
 	  {
          mutatedcodon = one + two + "T"; 
          dnainput2 = dnainput;
          dnanew = dnainput2.replace(i,3,mutatedcodon);
  	  dnalist.push_back(dnanew);
          dnanew = "";
          mutatedcodon = "";
          
          mutatedcodon = one + two + "C";           dnainput2 = dnainput;
          dnanew = dnainput2.replace(i,3,mutatedcodon);
  	  dnalist.push_back(dnanew);
          dnanew = "";
          mutatedcodon = "";
  	  
          mutatedcodon = one + two + "G"; 
          dnainput2 = dnainput;
          dnanew = dnainput2.replace(i,3,mutatedcodon);
  	  dnalist.push_back(dnanew);
          dnanew = "";          mutatedcodon = "";
    	  }

	if(cod[2] == 'T')
 	  {
          mutatedcodon = one + two + "A"; 
          dnainput2 = dnainput;
          dnanew = dnainput2.replace(i,3,mutatedcodon);
  	  dnalist.push_back(dnanew);
          dnanew = "";
          mutatedcodon = "";
          
          mutatedcodon = one + two + "C"; 
          dnainput2 = dnainput;
          dnanew = dnainput2.replace(i,3,mutatedcodon);
  	  dnalist.push_back(dnanew);
          dnanew = "";
          mutatedcodon = "";
  	  
          mutatedcodon = one + two + "G"; 
          dnainput2 = dnainput;
          dnanew = dnainput2.replace(i,3,mutatedcodon);
  	  dnalist.push_back(dnanew);
          dnanew = "";          mutatedcodon = "";
    	  }

    	if(cod[2] == 'G')
 	  {
          mutatedcodon = one + two + "A"; 
          dnainput2 = dnainput;
          dnanew = dnainput2.replace(i,3,mutatedcodon);
  	  dnalist.push_back(dnanew);
          dnanew = "";
          mutatedcodon = "";
          
          mutatedcodon = one + two + "C"; 
          dnainput2 = dnainput;
          dnanew = dnainput2.replace(i,3,mutatedcodon);
  	  dnalist.push_back(dnanew);
          dnanew = "";
          mutatedcodon = "";
  	  
          mutatedcodon = one + two + "T"; 
          dnainput2 = dnainput;
          dnanew = dnainput2.replace(i,3,mutatedcodon);
  	  dnalist.push_back(dnanew);
          dnanew = "";
          mutatedcodon = "";
    	  }

   	if(cod[2] == 'C')
 	  {
          mutatedcodon = one + two + "A"; 
          dnainput2 = dnainput;
          dnanew = dnainput2.replace(i,3,mutatedcodon);
  	  dnalist.push_back(dnanew);
          dnanew = "";
          mutatedcodon = "";
          
          mutatedcodon = one + two + "G"; 
          dnainput2 = dnainput;
          dnanew = dnainput2.replace(i,3,mutatedcodon);
  	  dnalist.push_back(dnanew);
          dnanew = "";
          mutatedcodon = "";
  	  
          mutatedcodon = one + two + "T"; 
          dnainput2 = dnainput;
          dnanew = dnainput2.replace(i,3,mutatedcodon);
  	  dnalist.push_back(dnanew);
          dnanew = "";
          mutatedcodon = "";
    	  }
    }
//cout<<"dnalist size = "<<dnalist.size()<<endl;
std::vector<string> mutprotein;std::vector<float> energyarray;
float energyoutput;

//convert dnalist into mutprotein list////////////////
for(int y=0; y<dnalist.size();y++)
   {	
    dnaseq = "";protein = "";
    dnaseq = dnalist[y];
    for(int z=0; z<(dnaseq.length()-2);z += 3)
       {
	codon = dnaseq.substr(z,3);
       	if(codon == "TAA"|| codon == "TAG" ||codon =="TGA")
	  {
	  break;
	  }
        else protein += translate(codon); 
   	}
    mutprotein.push_back(protein);
  }

    string proteinseq;
    string str;

    //cout<<"Number of mutated proteins = "<<mutprotein.size()<<endl;


//////////calculate an energy for each mutant protein///////////////////////////
    for(int h=0; h<mutprotein.size();h++)
        {
                proteinseq = "";
                proteinseq = mutprotein[h];
               

                if(proteinseq.length() < (dnainput.length()/3))
                  {
                  energyoutput = 0.0; 
                  //cout<<"Truncated protein "<<energyoutput<<endl;
                  outfile<<"Truncated protein "<<energyoutput<<endl;
                  //cout<<proteinseq<<endl;
		 //energyarray.push_back(energyoutput);
                  }  

		else if(proteinseq == protein)
                {
                energyoutput = origenergy; 
		//cout<<"Synonymous mutation"<<energyoutput<<endl;
                //outfile<<"Synonymous mutation: "<<endl;
                outfile<<proteinseq<<endl;
                outfile<<energyoutput<<endl;
                energyarray.push_back(energyoutput);
                }      

                else
                   {
                   ofstream outseq("mutatedprotein2");
                    outseq<<proteinseq<<endl;
                   //system("rm nurefined.pdb");

                   // This is the old code:: system("./scwrl3_lin/scwrl3 -i refined.pdb -o nurefined.pdb -s mutatedprotein > logfile");
                   //system("./scwrl3_lin/scwrl3 -i 1EGL.pdb -o nurefined2.pdb -s mutatedprotein2 > logfile");
                   
                   system("./scwrl4_lin/Scwrl4 -i 1EGL.pdb -o nurefined2.pdb -s mutatedprotein2 -h -v > logfile");
                   outseq.clear();
                   outseq.close();
                
                   remove("mutatedprotein2");
                   X.clear();
                   Y.clear();
                   Z.clear();
                   num.clear();
		   aaseq.clear();

                   line = "";
                   energyoutput = 0.0;

		   ifstream pdbfile;
		   pdbfile.open("nurefined2.pdb");                

                   if(!pdbfile)
                     {
                      cerr << "Unable to open pdb file";
                      exit(1);
                     }


 		   if(pdbfile.is_open())
  		     { 
                     //outfile<<"File is open"<<endl;
    		     while(! pdbfile.eof())
    		          {
      		          getline (pdbfile,line);
      		          int Epos = line.find("ATOM",0);
       		          int Epos2 = line.find("A",21);

       			   if(Epos2 != string::npos && Epos != string::npos)
       			      {
           	              // cout<<line<<endl;
                              str = "";
       			      str = line.substr(17,3);
      
     			      if(str == "ALA"){/*cout<<"A"<<endl;*/ aaseq.push_back('A');}
    			      if(str == "ARG"){/*cout<<"R"<<endl;*/ aaseq.push_back('R');}
       			      if(str == "ASN"){/*cout<<"N"<<endl;*/ aaseq.push_back('N');}
       		 	      if(str == "ASP"){/*cout<<"D"<<endl;*/ aaseq.push_back('D');}
       		 	      if(str == "CYS"){/*cout<<"C"<<endl;*/ aaseq.push_back('C');}
       		 	      if(str == "GLN"){/*cout<<"Q"<<endl;*/ aaseq.push_back('Q');}
      		 	      if(str == "GLU"){/*cout<<"E"<<endl;*/ aaseq.push_back('E');}
                 	      if(str == "GLY"){/*cout<<"G"<<endl;*/ aaseq.push_back('G');}
       			      if(str == "HIS"){/*cout<<"H"<<endl;*/ aaseq.push_back('H');}
       			      if(str == "ILE"){/*cout<<"I"<<endl;*/ aaseq.push_back('I');}
       			      if(str == "LEU"){/*cout<<"L"<<endl;*/ aaseq.push_back('L');}
      		       	      if(str == "LYS"){/*cout<<"K"<<endl;*/ aaseq.push_back('K');}
       		 	      if(str == "MET"){/*cout<<"M"<<endl;*/ aaseq.push_back('M');}
     		 	      if(str == "PHE"){/*cout<<"F"<<endl;*/ aaseq.push_back('F');}
      		  	      if(str == "PRO"){/*cout<<"P"<<endl;*/ aaseq.push_back('P');}
      			      if(str == "SER"){/*cout<<"S"<<endl;*/ aaseq.push_back('S');}
       			      if(str == "THR"){/*cout<<"T"<<endl;*/ aaseq.push_back('T');}
       			      if(str == "TRP"){/*cout<<"W"<<endl;*/ aaseq.push_back('W');}
       			      if(str == "TYR"){/*cout<<"Y"<<endl;*/ aaseq.push_back('Y');}
       			      if(str == "VAL"){/*cout<<"V"<<endl;*/ aaseq.push_back('V');}

 	       		      string seqnum = line.substr(23, 3);
   	      		      for(int i = 0; i < seqnum.length(); i++)
     			         {
    	      		         if(seqnum[i] == ' ')
        		            {
         		            seqnum.replace(i,1,"");
         		            i = i - 1;
       	    		            } 
       			         }
      			      int seqnumint = atoi(seqnum.c_str());
     			      num.push_back(seqnumint);

      	 		      string x = line.substr(31, 7);
         		      for(int i = 0; i < x.length(); i++)
        	                 {
          	                 if(x[i] == ' ')
          	                   {
          	                   x.replace(i,1,"");
          	                   i = i - 1;
          	                   } 
       	                         }
      			     float xint = atof(x.c_str());
      			     X.push_back(xint);

 
  	    	            string y = line.substr(39,7);
   	    	            for(int i = 0; i < y.length(); i++)
   	                       {
      	                        if(y[i] == ' ')
      	                          {
       	                           y.replace(i,1,"");
        	                   i = i - 1;
         	                  } 
                                }
    	       		    float yint = atof(y.c_str());
              	            Y.push_back(yint);


      			    string z = line.substr(47,7);
                  	    for(int i = 0; i < z.length(); i++)
        		       {
        	               if(z[i] == ' ')
          	                 {
         	                  z.replace(i,1,"");
          	                  i = i - 1;
          	                 } 
                               }
                           float zint = atof(z.c_str());
                           Z.push_back(zint);
          	           }
      	                 }
                   pdbfile.clear();
    	           pdbfile.close();
                   remove("nurefined2.pdb");
   	           }
             energyoutput = energy(X,Y,Z,num,aaseq.size(),aaseq, contactenergy);
             energyarray.push_back(energyoutput);
	     perlinput << energyoutput << endl;
             }
      }

//cout<<"Size of energy array = "<<energyarray.size()<<endl;

float ddG = 0.0;
float avddG;

for(int b=0; b<energyarray.size(); b++)   {    ddG = ddG + abs(origenergy - energyarray[b]);     }


avddG = ddG / energyarray.size();
robustout<<avddG<<endl; 
perlinput<<"robustness ="<<avddG<<endl;

system("perl countstable.pl perlinput >> ./stabdestab");
system("rm perlinput");

}


///////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////this is the main procedure/////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////

int main(int argc, char *argv[])
{


/////////the following refines the initial pdb file////////////////////////////////////////////
remove("refined.pdb");
remove("stabdestab");
//system("./scwrl3_lin/scwrl3 -i protonlynu.pdb -o refined.pdb -s aa > logfile");
//system("./scwrl3_lin/scwrl3 -i 1EGL.pdb -o refined.pdb -s aa > logfile");




system("./scwrl4_lin/Scwrl4 -i 1EGL.pdb -o refined.pdb -s aa -h -v > logfile1");

//////the following extracts information from the pdb file ie. x,y,z, amino acid seq and numbers///

ifstream pdbfile;
pdbfile.open("refined.pdb");
std::vector<int> num;
std::vector<float> X;
std::vector<float> Y;
std::vector<float> Z;
std::vector<char> aaseq;
string line;

  if (pdbfile.is_open())
  {
    while (! pdbfile.eof())
    {
      getline (pdbfile,line);
      int Epos = line.find("ATOM",0);
      int Epos2 = line.find("A",21);

      
       if(Epos2 != string::npos && Epos != string::npos)
       {
      // cout<<line<<endl;
       string str = line.substr(17,3);
      
       if(str == "ALA"){/*cout<<"A"<<endl;*/ aaseq.push_back('A');}
       if(str == "ARG"){/*cout<<"R"<<endl;*/ aaseq.push_back('R');}
       if(str == "ASN"){/*cout<<"N"<<endl;*/ aaseq.push_back('N');}
       if(str == "ASP"){/*cout<<"D"<<endl;*/ aaseq.push_back('D');}
       if(str == "CYS"){/*cout<<"C"<<endl;*/ aaseq.push_back('C');}
       if(str == "GLN"){/*cout<<"Q"<<endl;*/ aaseq.push_back('Q');}
       if(str == "GLU"){/*cout<<"E"<<endl;*/ aaseq.push_back('E');}
       if(str == "GLY"){/*cout<<"G"<<endl;*/ aaseq.push_back('G');}
       if(str == "HIS"){/*cout<<"H"<<endl;*/ aaseq.push_back('H');}
       if(str == "ILE"){/*cout<<"I"<<endl;*/ aaseq.push_back('I');}
       if(str == "LEU"){/*cout<<"L"<<endl;*/ aaseq.push_back('L');}
       if(str == "LYS"){/*cout<<"K"<<endl;*/ aaseq.push_back('K');}
       if(str == "MET"){/*cout<<"M"<<endl;*/ aaseq.push_back('M');}
       if(str == "PHE"){/*cout<<"F"<<endl;*/ aaseq.push_back('F');}
       if(str == "PRO"){/*cout<<"P"<<endl;*/ aaseq.push_back('P');}
       if(str == "SER"){/*cout<<"S"<<endl;*/ aaseq.push_back('S');}
       if(str == "THR"){/*cout<<"T"<<endl;*/ aaseq.push_back('T');}
       if(str == "TRP"){/*cout<<"W"<<endl;*/ aaseq.push_back('W');}
       if(str == "TYR"){/*cout<<"Y"<<endl;*/ aaseq.push_back('Y');}
       if(str == "VAL"){/*cout<<"V"<<endl;*/ aaseq.push_back('V');}

 	       string seqnum = line.substr(23, 3);
   	       for (int i = 0; i < seqnum.length(); i++)
     		     {
    	          if (seqnum[i] == ' ')
        	      {
         	     seqnum.replace(i,1,"");
         	     i = i - 1;
       	    		   } 
       		     }
      		int seqnumint = atoi(seqnum.c_str());
     		num.push_back(seqnumint);


      	 	string x = line.substr(31, 7);
         	 for (int i = 0; i < x.length(); i++)
        	  {
          	    if (x[i] == ' ')
          	    {
          	    x.replace(i,1,"");
          	    i = i - 1;
          	    } 
       	          }
      		float xint = atof(x.c_str());
      		X.push_back(xint);

 
  	       string y = line.substr(39,7);
   	       for (int i = 0; i < y.length(); i++)
   	       {
      	        if (y[i] == ' ')
      	        {
       	         y.replace(i,1,"");
        	 i = i - 1;
         	 } 
               }
    	       float yint = atof(y.c_str());
               Y.push_back(yint);


      		string z = line.substr(47,7);
         	 for (int i = 0; i < z.length(); i++)
        	  {
        	     if (z[i] == ' ')
          	     {
         	     z.replace(i,1,"");
          	     i = i - 1;
          	     } 
                  }
                float zint = atof(z.c_str());
                Z.push_back(zint);

               // cout<<seqnum<<"x"<<endl;
          }
      }
    pdbfile.close();

   }

/////////////////////the following calls the energy calculation///////////

ifstream energyfile;
energyfile.open("vendruscolomatrix"); 
float contactenergy[20][20]; 

   	for (int i=0;i<20;i++)
    	{
		for(int j=0;j<20;j++)
    		{
      			energyfile>>contactenergy[i][j];
    		}
	}

energyfile.close();

float origenergyoutput;

origenergyoutput = energy(X,Y,Z,num,aaseq.size(),aaseq, contactenergy);

outfile<<"initial refined energy = "<<origenergyoutput<<endl;
robustout<<"initial energy = "<<origenergyoutput<<endl;


///////////////////////////////the next part is the evolution part/////////////////////////////

ifstream dnafile;
dnafile.open("dna");
string dnainput, codon, out, mutatedprot, stableprot;
dnafile>>dnainput;
outfile<<"DNA Seq : "<<dnainput<<endl;
dnafile.close();
float energyoutput;

std::vector<string> dna;
std::vector<string> mutprotein;
std::vector<string> stable;
//std::vector<string> population;

int m[POPN]; int z;
StochasticLib1 stochastic(seed + rg.IRandom(0,1000)); 

for(int i=0;i<POPN;i++)
{
	dna.push_back(dnainput);
}

for(int a=0;a<GEN;a++)
    {
     
        for(int i=0;i<POPN;i++)
           {
/////////////////////the following generates the mutated protein////////////////////////////////
		//m[i] = stochastic.Poisson(0.14); /////this is the number of mutations per sequence and gives 0.002 per nucleotide per generation /////////////
                m[i] = stochastic.Poisson(1); 
              //  m[i] = stochastic.Poisson(0.1); 
                out = "";
                out = Mutate(dna[i], m[i]); 
                outfile<<"Number of mutations = "<<m[i]<<endl;

                //cout<<out<<endl;
                mutatedprot = "";
		for(int j=0;j<out.length()-2;j += 3)
		{
			codon = out.substr(j,3);
       			if(codon == "TAA"|| codon == "TAG" ||codon =="TGA")
			{
				break;
			}

			else mutatedprot += translate(codon);
		}
               // cout<<mutatedprot<<endl;
                outfile<<mutatedprot<<endl;
                outfile<<i<<endl;

///////////////the following calculates the energy of the mutated protein///////////////////////

                int length = mutatedprot.length();
                if(length < 70)
                {
                energyoutput = 0.0;
                }  

                else
                {
                //remove("nurefined.pdb");
                ofstream outseq("mutatedprotein");
                outseq<<mutatedprot<<endl;
                //system("./scwrl3_lin/scwrl3 -i 1EGL.pdb -o nurefined.pdb -s mutatedprotein > logfile");
//if(i==0){system("./scwrl4_lin/Scwrl4 -i 1EGL.pdb -o nurefinedtest.pdb -s mutatedprotein -h > logfile2");}
                system("./scwrl4_lin/Scwrl4 -i 1EGL.pdb -o nurefined.pdb -s mutatedprotein -h -v  > logfile");
                outseq.close();             
                remove("mutatedprotein");

                X.clear();
                Y.clear();
                Z.clear();
                num.clear();
		aaseq.clear();
                line = "";
                energyoutput = 0.0;
            
		ifstream pdbfile;
		pdbfile.open("nurefined.pdb");

 		 if(pdbfile.is_open())
  		   {
    		  while(! pdbfile.eof())
    		     {
      		    getline (pdbfile,line);
      		    int Epos = line.find("ATOM",0);
       		    int Epos2 = line.find("A",21);

      
       			if(Epos2 != string::npos && Epos != string::npos)
       			  {
           	          // cout<<line<<endl;
       			  string str = line.substr(17,3);
      
     			  if(str == "ALA"){/*cout<<"A"<<endl;*/ aaseq.push_back('A');}
    			  if(str == "ARG"){/*cout<<"R"<<endl;*/ aaseq.push_back('R');}
       			  if(str == "ASN"){/*cout<<"N"<<endl;*/ aaseq.push_back('N');}
       		 	  if(str == "ASP"){/*cout<<"D"<<endl;*/ aaseq.push_back('D');}
       		 	  if(str == "CYS"){/*cout<<"C"<<endl;*/ aaseq.push_back('C');}
       		 	  if(str == "GLN"){/*cout<<"Q"<<endl;*/ aaseq.push_back('Q');}
      		 	  if(str == "GLU"){/*cout<<"E"<<endl;*/ aaseq.push_back('E');}
                 	  if(str == "GLY"){/*cout<<"G"<<endl;*/ aaseq.push_back('G');}
       			  if(str == "HIS"){/*cout<<"H"<<endl;*/ aaseq.push_back('H');}
       			  if(str == "ILE"){/*cout<<"I"<<endl;*/ aaseq.push_back('I');}
       			  if(str == "LEU"){/*cout<<"L"<<endl;*/ aaseq.push_back('L');}
      			  if(str == "LYS"){/*cout<<"K"<<endl;*/ aaseq.push_back('K');}
       		 	  if(str == "MET"){/*cout<<"M"<<endl;*/ aaseq.push_back('M');}
      		 	  if(str == "PHE"){/*cout<<"F"<<endl;*/ aaseq.push_back('F');}
      		  	  if(str == "PRO"){/*cout<<"P"<<endl;*/ aaseq.push_back('P');}
      			  if(str == "SER"){/*cout<<"S"<<endl;*/ aaseq.push_back('S');}
       			  if(str == "THR"){/*cout<<"T"<<endl;*/ aaseq.push_back('T');}
       			  if(str == "TRP"){/*cout<<"W"<<endl;*/ aaseq.push_back('W');}
       			  if(str == "TYR"){/*cout<<"Y"<<endl;*/ aaseq.push_back('Y');}
       			  if(str == "VAL"){/*cout<<"V"<<endl;*/ aaseq.push_back('V');}

 	       		string seqnum = line.substr(23, 3);
   	      		 for(int i = 0; i < seqnum.length(); i++)
     			    {
    	      		    if(seqnum[i] == ' ')
        		      {
         		      seqnum.replace(i,1,"");
         		      i = i - 1;
       	    		      } 
       			     }
      			int seqnumint = atoi(seqnum.c_str());
     			num.push_back(seqnumint);


      	 		string x = line.substr(31, 7);
         		for(int i = 0; i < x.length(); i++)
        	           {
          	          if(x[i] == ' ')
          	            {
          	            x.replace(i,1,"");
          	            i = i - 1;
          	            } 
       	                   }
      			float xint = atof(x.c_str());
      			X.push_back(xint);

 
  	    	        string y = line.substr(39,7);
   	    	        for(int i = 0; i < y.length(); i++)
   	                   {
      	                   if(y[i] == ' ')
      	                     {
       	                     y.replace(i,1,"");
        	             i = i - 1;
         	             } 
                           }
    	       		float yint = atof(y.c_str());
              	        Y.push_back(yint);


      			string z = line.substr(47,7);
                  	 for(int i = 0; i < z.length(); i++)
        		    {
        	            if(z[i] == ' ')
          	              {
         	              z.replace(i,1,"");
          	              i = i - 1;
          	              } 
                            }
                        float zint = atof(z.c_str());
                        Z.push_back(zint);

                        // cout<<seqnum<<"x"<<endl;
          	}
      	    }
    	pdbfile.close();
        remove("nurefined.pdb"); // removing this line speeds up the code considerably
   	}

        energyoutput = energy(X,Y,Z,num, aaseq.size(), aaseq, contactenergy);
        
        }
        outfile<<"energy = "<<energyoutput<<endl;

        
//////////the following is the population drift and selection model/////////////////////////
              // if(energyoutput==origenergyoutput)
	     //	if(energyoutput<=origenergyoutput+0235) 
               if(energyoutput < 0)
            //   if((energyoutput<=origenergyoutput+0.31) && (energyoutput >=origenergyoutput-0.31))

               {
		stable.push_back(out);
                outfile<<"this is stablesize "<<stable.size()<<endl;
		}
        else{}
        
    }


/////the population needs to be randomly re-populated so that the size is equivalent to POPN////
       
        int initialsize = stable.size();
        outfile<<"\nGeneration = "<<a + 1<<endl;
        outfile<<"Selected popn = "<<stable.size()<<endl;
           
        if(stable.size() != 0)
          {
           for(int x=0; x<POPN-initialsize; x++)
              {int z = rg.IRandom(0,initialsize-1);
              stable.push_back(stable[z]);
              //outfile<<"stable[z] = "<<stable[z]<<endl;
              }
 
           outfile<<"Final popn = "<<stable.size()<<endl;
        ////transfer stable to a new vector, 'dna'
           dna.clear();

           outfile<<"\nFINAL POPULATION (AFTER SELECTION)\n"<<endl;
           ofstream popnout("popndnaseqs");
           for(int b=0; b<stable.size();b++)
              {dna.push_back(stable[b]);
              outfile<<stable[b]<<endl;   
              popnout<<stable[b]<<endl;
              }
           
           stable.clear();
 
          // stable.erase(stable.begin(), stable.end());
          // if(stable.empty())
          //   {outfile<<"EMPTY!"<<endl;}
          //    outfile<<"stable zero = "<<stable[9]<<"\n"<<endl;
        

//here the stable protein sequences in the selected population are printed/////////////////////

          for(int c=0; c<dna.size();c++)
             {
              stableprot = "";
              for(int j=0;j<dna[c].length()-2;j += 3)
		 {
		 codon = dna[c].substr(j,3);
       	  	 stableprot += translate(codon);       
		 }
              outfile<<stableprot<<endl;
              }
          outfile<<"\n"<<endl; 
          }
     
          else{
               for(int c=0; c<dna.size();c++)
                  {
                   stableprot = "";
                   for(int j=0;j<dna[c].length()-2;j += 3)
		       {
		       codon = dna[c].substr(j,3);
       	               stableprot += translate(codon);   
		       }
                   outfile<<stableprot<<endl;
                   }
              outfile<<"\n"<<endl;
              }

  robustness(contactenergy);
  remove("mostcommonout");
  remove("popndnaseqs");
  }

return 0;
}







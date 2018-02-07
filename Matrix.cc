//@author Andrew McKeighan
#include <iostream> // provides objects like cin and cout for sending data
                    // to and from the standard streams input and output.
		    // These objects are part of the std namespace.
#include <cstdlib>  // has exit etc.
#include <fstream>  // file streams and operations
#include <sstream>  // string streams and operations
using namespace std; // a container for a set of identifiers.
                     // Because of this line, there is no need to use std::cout
#include <typeinfo> // This header file is included for using typeid.
#include <stdexcept>

#include <stdio.h>
#include <string.h>
#include <ctype.h>  // A library with functions for testing characters.


#include "wrapper.hh"
#include "Alignment.hh"
#include "Direct.hh"
#include "Matrix.hh"
#include "Encoded.hh"
#include "prototype.hh"


Matrix:: Matrix(Direct &seqone, Direct &seqtwo, struct scoretp &param)
// Input parameter seqone is a reference to a Direct object holding one DNA sequence.
// Input parameter seqtwo is a reference to a Direct object holding another DNA sequence.
// Input parameter param is a reference to a struct scoretp variable holding the scoring parameters.

// The constructor sets the member reference origin to seqone, and
// the member reference derived to seqtwo,
// through the initialization list before the constructor body.
//
// The constructor sets the member rowind to the return value from getLength() of seqone.
// The constructor sets the member colind to the return value from getLength() of seqtwo.
// If both rowind and colind are zero, then prints an error message.
// Then it allocates space for each
// matrix in the mattp structure by calling get2dspace(). Next it computes,
// in reverse order, the three matrices by the dynamic programming algorithm
// for aligning the two sequences pointed to by getLength() of seqone
// that of seqtwo. The sequence in seqone is treated as the sequence A.
: origin(seqone), derived(seqtwo)  // an initialization list.
{
  rowind = origin.getLength();
  colind = derived.getLength();


  if(rowind == 0 && colind ==0){
    fatal("rowind and colind both equal 0");
  }
  Dmat = get2dspace(rowind, colind);
  Imat = get2dspace(rowind, colind);
  Smat = get2dspace(rowind, colind);


  Smat[rowind][colind] = 0;
  Dmat[rowind][colind] = Smat[rowind][colind]-(param.gopen);
  Imat[rowind][colind] = Smat[rowind][colind]-(param.gopen);

  for(int j = (colind-1); j >= 0; --j){
      Imat[rowind][j] = Imat[rowind][j+1]-(param.gext);
      Smat[rowind][j] = Imat[rowind][j];
      Dmat[rowind][j] = Smat[rowind][j]-(param.gopen);
  }
  for(int i = (rowind-1); i >= 0; --i){
      Dmat[i][colind] = Dmat[i+1][colind]-(param.gext);
      Smat[i][colind] = Dmat[i][colind];
      Imat[i][colind] = Smat[i][colind]-(param.gopen);

      for(int j = (colind)-1; j >= 0; --j){
          Dmat[i][j] = max(Dmat[i+1][j]-(param.gext), Smat[i+1][j]-(param.gext)-(param.gopen));
          Imat[i][j] = max(Imat[i][j+1]-(param.gext), Smat[i][j+1]-(param.gext)-(param.gopen));
          Smat[i][j] = trimax(Smat[i+1][j+1] + compare(origin.getSeq()[i+1], derived.getSeq()[j+1], param), Dmat[i][j], Imat[i][j]);
    }
}


}
int Matrix:: max(int one, int two){
    if(one >= two){
        return one;
    }else{
        return two;
    }
}
int Matrix:: compare(char a, char b, struct scoretp &param){
    if(a == b){
        return param.match;
    }else{
        return param.mismat;
    }
}
int Matrix:: trimax(int one, int two, int three){
    return max(max(one, two) , three);
}

//}

// Frees heap memory.
Matrix:: ~Matrix()
{
  free2dspace(rowind, Dmat);
  free2dspace(rowind, Imat);
  free2dspace(rowind, Smat);
}

int** Matrix:: get2dspace(int rowind, int colind)
// Allocates a 2-dimensional int array of (rowind + 1) by (colind + 1)
// and returns its memory address.
// Prints out an error message if rowind or colind is less than 0.
{
  int** arr = new int*[rowind + 1];
  for(int i = 0; i < (rowind+1); ++i){
    arr[i] = new int[colind + 1];
  }
  return arr;
}

void Matrix:: free2dspace(int rowind, int **arr)
// Deallocates a 2-dimensional int array of first dimension 
// (rowind + 1) by releasing the memory for arr[j] for each j
// and then releasing the memory for arr.
// Prints out an error message if rowind is less than 0.
{
  if(rowind < 0){
    fatal("rowind < 0");
  }
  for(int i = 0; i < (rowind + 1); ++i){
    delete[] arr[i];
    arr[i]=NULL;
  }
  delete[] arr;
      
}

Direct& Matrix:: getOrigin() const
// Returns origin.
{
  return origin;
}

Direct& Matrix:: getDerived() const
// Returns derived.
{
  return derived;
}

int Matrix:: getRowInd() const
// Returns rowind.
{
  return rowind;
}

int Matrix:: getColInd() const
// Returns colind.
{
  return colind;
}

int** Matrix:: getMat(char kind) const
// If kind is 'S', then returns Smat.
// If kind is 'D', then returns Dmat.
// If kind is 'I', then returns Imat.
// Otherwise, prints out an error message.
{
  if(kind == 'S'){
    return Smat;
  }else if(kind == 'I'){
    return Imat;
  }else if(kind == 'D'){
    return Dmat;
  }else{
    fatal("getMat did not recieve a proper kind");
  }
  return 0;
}

string Matrix:: toString(char kind) const
// Input parameter kind denotes matrix type: 'D', 'I', or 'S'.
// The function returns the matrix type and each value in the matrix
// in the form of a string.
// This function is for your own use to check on each matrix,
// so any format is OK.
{
  int i = 0;
  int j = 0;
  ostringstream mato;
  

    if(kind == 'S' || kind == 'D' || kind == 'I'){
        for(i = 0; i<=origin.getLength(); ++i){
            for(j=0; j<=derived.getLength(); ++j){
              mato.width(4);
              mato << getMat(kind)[i][j];
                //cout<< getMat(kind)[i][j] <<" ";
            }
            mato << "\n";
            //printf("\n");
        }
    }
    return mato.str();

}

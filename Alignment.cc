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
#include "prototype.hh"


Alignment:: Alignment(Matrix &matobj, struct scoretp &param)
// Input parameter matobj is a reference to a Matrix object.
// Input parameter param is a reference to a struct scoretp variable with scoring parameters.

// The constructor initializes some member variables by calling public
// functions of Matrix from the object matobj.
// Then it allocates memory for each of the three char arrays,
// produces an optimal alignment by tracing through the matrices in the object matobj,
// and saves the alignment along with its score and length and editnum
// and subinsertlen in corresponding member variables.
//
// The alignment is represented by using three rows of characters (three char arrays).
// If the lengths of both sequences is zero, prints out an error message.

// Note that the length of an optimal alignment cannot exceed the sum of
// the lengths of the two DNA sequences.
: origin(matobj.getOrigin()), derived(matobj.getDerived()){
  int i = 0;
  int j = 0;
  int k = 0;
  int** S = matobj.getMat('S');
  int** D = matobj.getMat('D');
  int** I = matobj.getMat('I');
  alignlen = '0';
  //int** mat = S;
  char matc = 'S';
  int m = matobj.getRowInd();
  int n = matobj.getColInd();
  int intlen = 0;
  top = new char[m + n + 1];
  mid = new char[m + n + 1];
  bot = new char[m + n + 1];
  for(int y = 0; y < m+n+1; ++y){
    top[y] = ' ';
    mid[y] = ' ';
    bot[y] = ' ';
  }



  int q = param.gopen;
  int r = param.gext;
  editnum = 0;
  subinsertlen = 0;
  intlen = 0;

  while(i <= m || j <= n){
    if(matc == 'S'){
      if(i == m && j == n){
        break;
      }
      if(j == n || S[i][j] == D[i][j]){
        //mat = D;
        matc = 'D';
        editnum++;
        continue;
      }
      if(i == m || S[i][j] == I[i][j]){
        //mat = I;
        matc = 'I';
        editnum++;
        intlen = 0;
        continue;
      }
      top[k] = origin.getSeq()[i+1];
      bot[k] = derived.getSeq()[j+1];
      if(origin.getSeq()[i+1] == derived.getSeq()[j+1]){
        mid[k] = '|';
      }else{
        mid[k] = ' ';
      }

      i++;
      j++;
      k++;
      if(origin.getSeq()[i] != derived.getSeq()[j]){ // do i need getOrigin()
        editnum++;
        subinsertlen++;
      }
      continue;
    }
    if(matc == 'D'){
      top[k] = origin.getSeq()[i+1];
      mid[k] = '-';
      bot[k] = ' ';

      if(i == (m-1) || D[i][j] == (S[i+1][j] - q - r)){
        //mat = S;
        matc = 'S';
      }
      i++;
      k++;
      continue;
    }
    if(matc == 'I'){
      top[k] = ' ';
      mid[k] = '-';
      bot[k] = derived.getSeq()[j+1];
      intlen++;
      if(j == (n-1) || I[i][j] == (S[i][j+1] - q - r)){
        //mat = S;
        matc = 'S';
        subinsertlen += intlen;
      }
      ++k;
      ++j;
      continue;
    }
  }
alignlen = k;
top[k+1] = '\0';//it was k+1
bot[k+1] = '\0';
mid[k+1] = '\0';
score = S[0][0];

}

// Frees heap memory.
Alignment::  ~Alignment() // destructor definition
{
  delete[] top;
  delete[] mid;
  delete[] bot;
  top = NULL;
  mid = NULL;
  bot = NULL;
}

Direct& Alignment:: getOrigin() const
// returns origin.
{
  return origin;
}

Direct& Alignment:: getDerived() const
// returns derived.
{
  return derived;
}

int Alignment:: getScore() const
// returns alignment score.
{
  return score;
}

int Alignment:: getEditNum() const
// returns editnum.
{
  return editnum;
}

int Alignment:: getSubInsertLen() const
// Returns subinsertlen.
{
  return subinsertlen;
}

int Alignment:: getAlignLen() const
// Returns alignment length.
{
  return alignlen;
}

char* Alignment:: getRow(char kind) const
// Returns top if kind is 'T'.
// Returns mid if kind is 'M'.
// Returns bot if kind is 'B'.
// Prints out an error if kind has no expeced value.
{
  if(kind == 'T'){
    return top;
  }else if(kind == 'M'){
    return mid;
  }else if(kind == 'B'){
    return bot;
  }else{
    fatal("getRow has an improper value");
  }
  return 0;
}

string Alignment:: toString() const
// Input parameter sonept is a pointer to a structure with the features of one DNA sequence.
// Input parameter stwopt is for another DNA sequence.
// Input parameter matspt is a pointer to a structure with three matrices.
// Input parameter algopt is a pointer to a structure for an optimal alignment.

// The function returns a summary of sequence and alignment information
// and the alignment in the form of string.
// The summary includes the name and length of each sequence,
// the score and length of the alignment. the number of edit operations,
// and the length of substitutions and insertions.
// The alignment is reported in sections of 70 characters, with each section consisting
// of three rows. The sequence positions of the first DNA letters in each section are
// reported in the left margin of 10 spaces.
{
  ostringstream fin;
  ostringstream topln;
  ostringstream fubar;
  ostringstream botln;
  int i = 0;
  int j = 0;
  int row = (alignlen/70)+1;

  for(i=0; i< row; ++i){
        for(j=0; j<70; ++j){
            if(mid[j] != '\0'){
              topln << top[j];
              fubar << mid[j];
              botln << bot[j];
            }else{
              break;
            }
        }
        fin << topln.str() << '\n';
        fin << fubar.str() << '\n';
        fin << botln.str() << '\n';
    }

  return fin.str();
}

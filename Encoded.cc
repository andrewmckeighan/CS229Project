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

Encoded:: Encoded(Alignment &obj)
// Input parameter obj is a reference to an Alignment object holding an optimal
// alignment along with sequence and other information.

// The constructor initializes some members by calling public functions
// of Alignment from the object obj. Among these members are
// origin, the reference to the original sequence,
// editnum, the size of the array operation, and
// subinsertlen, the size of the array subinsertion.
// the name and length of the derived sequence.
// The constructor allocates memory to operation and subinsertion.
// Then it constructs the two arrays operation and subinsertion
// by follwing the method in Figure 2 of the project description.
// The construction uses top & mid & bot, the three rows of the optimal alignment,
// which are provided by the public functions of Alignment from the object obj.
: origin(obj.getOrigin()) {

  subinsertlen = obj.getSubInsertLen();
  editnum = obj.getEditNum();
  int i = 1;
  int j = 1;
  int opind = 0;
  int snind = 0;
  int gaplen = 0;
  int ind = 0;
  dname = obj.getDerived().getName();
  dlength = obj.getDerived().getLength();
  operation = new struct Edit[editnum];
  subinsertion = new char[subinsertlen+1];

  for (ind = 0; ind < obj.getAlignLen();){
    if(obj.getRow('M')[ind] == '|'){
      i++;
      j++;
      ind++;
      continue;
    }
    if(obj.getRow('M')[ind] == ' '){
      operation[opind].position = i;
      operation[opind++].indel = 0;
      subinsertion[snind++] = obj.getRow('B')[ind];
      i++;
      j++;
      ind++;
      continue;
    }
    int k = ind;
    while(obj.getRow('M')[k] == '-'){
      gaplen++;
      k++;
    }

    operation[opind].position = i;
    if(obj.getRow('M')[ind] == '-' && obj.getRow('B')[ind] == ' '){
      i += gaplen;
      operation[opind++].indel = -gaplen;
      ind += gaplen;
    }else{
      j += gaplen;
      operation[opind++].indel = gaplen;
      for( ; gaplen; gaplen--){
        subinsertion[snind++] = obj.getRow('B')[ind++];
      }//insertion gap
    }
  }
  


}

Encoded:: ~Encoded()
// Frees heap memory.
{
  delete[] operation;
  delete[] subinsertion;
  operation = NULL;
  subinsertion = NULL;
}

char* Encoded:: getDSeq() const
// Deirves a sequence of length dlength and returns it in a char array
// of length dlength + 2. The sequence is in the array from index 1
// to index dlength, with the element at index 0 set to ' '
// and that at index dlength + 1 to '\0'.
// The function allocates memory to the array and generates
// the sequence B by using the method in Figure 3 of the project description.
// The information on the sequence A is obtained by calling the public
// functions of Direct from the object reference origin.
{
  int i = 1;
  int j = 1;
  char* A = origin.getSeq();
  char* B = new char[dlength+1];
  B[0] = ' ';
  int snind = 0;
  int opind = 0;
  int oppos = 0;
  short indel = 0;
  int k = 0;

  

  for( ; opind < editnum; opind++){
    oppos = operation[opind].position;
    indel = operation[opind].indel;
    if(i < oppos){
      for(k = 0; k <= (oppos-1-i); k++){ // make sure that it is <= or <
        B[j+k] = A[i+k];
        cout<< "B "<< B[j+k]<< " A "<<A[i+k]<<endl;
      }
      i = oppos;
      j += (oppos-i);
    }
    if(indel < 0){
      i -= indel;
    }else{
      if(! indel){
        i++;
      }
      B[j++] = subinsertion[snind++];
      for( ; indel > 1; indel--){
        B[j++] = subinsertion[snind++];
      }
    }
//    if(i <= origin.getLength()){
//      for(k = 0; k < origin.getLength()-i; k++){ // make sure that it is <= or <
//        B[j+k] = A[i+k];
//      }
//    }
//    B[dlength+1] = '\0';
  }

  if(i <= origin.getLength()){
      for(k = 0; k <= origin.getLength()-i; k++){ // make sure that it is <= or <
        B[j+k] = A[i+k];
      }
    }
    B[dlength+1] = '\0';
  cout<<B<<endl;
  return B;
}

Direct& Encoded:: getOrigin() const
// Returns origin.
{
 return origin;
}

int Encoded:: getEditNum() const
// Returns editnum.
{
  return editnum;
}

struct Edit* Encoded:: getOperation() const 
// Returns operation.
{
  return operation;
}

int Encoded:: getSubInsertLen() const
// Returns subinsertlen.
{
  return subinsertlen;
}

char* Encoded:: getSubInsertion() const
// Returns subinsertion.
{
  return subinsertion;
}

int Encoded:: getDLength() const
// Returns dlength.
{
  return dlength;
}

string Encoded:: getDName() const
// Returns dname.
{
  return dname;
}

string Encoded:: toString() const
// Generates and returns a string form of its contents.
// The information includes the following lines:
// Name of the encoded sequence: ...
// Length of the encoded sequence: ...
// Number of edit operations: ...
// Length of substitutions and insertions: ...
// Concatenation of subs and inserts in order: ...
//
// and each of the operations in the array operation.
{
string temp = "";
int chk = 0;

ostringstream fin;
fin<<"Name of the encoded sequence: "<< dname << "\n";
fin<<"Length of the encoded sequence: "<< dlength << "\n";
fin<<"Number of edit operations: "<< editnum << "\n";
fin<<"Length of substitutions and insertions: "<< subinsertlen << "\n";
fin<<"Concatenation of subs and inserts in order: " << "\n";

for(int i = 0; i <= subinsertlen; ++i){
  chk = 0;
  if(operation[i].indel > 0){
    temp = "Insertion: ";
  }else if(operation[i].indel == 0){
    temp = "Substitution: ";
  }else{
    temp = "Deletion";
    chk = 1;
  }
  fin<< "Operation "<< i << "::" << "  Position: "<< operation[i].position << "  Indel: "
  << operation[i].indel << "  " << temp;
  if(chk == 0){
     fin<< subinsertion[i]<< endl;
  }else{
    fin<<endl;
  }
}


return fin.str();
}

bool Encoded:: operator<=(Encoded &rightobj) const
{
  if(getNumDiff() < rightobj.getNumDiff()){
    return true;
  }else if(getNumDiff() > rightobj.getNumDiff()){
    return false;
  }else{
    if(subinsertlen <= rightobj.getSubInsertLen()){
      return true;
    }else{
      return false;
    }
  }
}

int Encoded:: getNumDiff() const
{
  int i = 0;
  int tlen = 0;
  for(i = 0; i < editnum; i++){
    if(operation[i].indel < 0){
      tlen += -operation[i].indel;
    }else if(operation[i].indel > 0){
      tlen += operation[i].indel;
    }else{
      tlen++;
    }
  }
  return tlen;
}

Compressed:: Compressed(Alignment &obj) : Encoded(obj)
{

}

bool Compressed:: operator<=(Encoded &rightobj) const
{
  int i = 0;
  if(getNumDiff() < rightobj.getNumDiff()){
    return true;
  }else if(getNumDiff() > rightobj.getNumDiff()){
    return false;
  }else{
    if(getSubInsertLen() < rightobj.getSubInsertLen()){
      return true;
    }else if(getSubInsertLen() > rightobj.getSubInsertLen()){
      return false;
    }else{
      while(i <= getSubInsertLen()){
        if(getSubInsertion()[i] != rightobj.getSubInsertion()[i]){
          return false;
        }
        i++;
      }
      return true;
    }
  }
}


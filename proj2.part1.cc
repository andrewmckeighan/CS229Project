// Template for Part 1 of Project 2
// @author Andrew McKeighan

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

void fatal(const char *msg); // a wrapper function
void fatal(const char *msg, const char *name); // a wrapper function
void ckopeninf(ifstream &infile, const char *fname); // a wrapper function

struct scoretp  // a structure for scoring parameters
 { int  match;  // a positive score for a pair of identical DNA letters   
   int  mismat; // a negative score for a pair of different DNA letters
   int  gopen;  // a negative score for a gap
   int  gext;   // a negative score for each letter in a gap
 };

struct Edit        // a structure for an edit operation
 { int   position; // a sequence position to which the edit operation is applied
   short indel;    // indel < 0, the edit is a deletion gap of length -indel;
                   // indel > 0, the edit is an insertion gap of length indel;
                   // indel = 0, the edit is a substitution.
 };

class Direct   // a class for an object holding a sequence
 { private:
    string name;   // the name of the sequence
    int    length; // the length of the sequence
    char   *seq;   // an array holding the sequence from index 1 to index length
                   // with seq[0] set to ' ' and seq[length+1] set to '\0'.
		   // The length of seq is length + 2.

   public:
    Direct(const char *fname); // normal constructor
    Direct(string &stag, int slen, char *sarr); // an alternative constructor
    Direct(const Direct &obj); // copy constructor
    Direct();                  // default constructor
    ~Direct();                 // destructor
    string getName() const;    // returns name.
    int    getLength() const;  // returns length.
    char*  getSeq() const;     // returns seq.
 };

Direct:: Direct(const char *fname)
// Input parameter fname is a pointer to a string for storing the name of a DNA sequence file.
// A Direct object is constructed to save the name,
// length and DNA letters of the sequence in the file.
// A proper amount of memory needs to be allocated for the seq member.

// This constructor opens the file whose name is pointed to by fname,
// reads the name, length and DNA letters of the sequence in the file,
// dynamically allocates a proper amount of memory for the seq member
// and saves the data into the member variables.
// The constructor prints out a proper error message if fname is NULL,
// no file exists with the given name.
{
  ifstream fp;
//  fp.open(fname);//*
  ckopeninf(fp, fname);
  char next;
  int chk = 0;
  int namelength = 0;
  int seqlength = 0;
  int i = 0;
  int j = 0;

  if(fname = NULL){
    fatal("fname is NULL.");
  }
  next = fp.get();
  while((next = fp.get()) != EOF && next != '\n'){
    ++namelength;
  }
  while((next = fp.get()) != EOF){
    ++seqlength;
  }
    //if(next == '\n'){
    //  chk = 1;
    //}
    //if(next != '\n' && chk == 0){
    //  ++namelength;
    //}
    //if(next != '\n' && chk == 0){
    //  ++seqlength;
    //}
    //next = fp.get();
  //}
  this->length = seqlength -1;

  seq = new char[seqlength + 2];
  fp.clear();
  fp.seekg(0);
  chk = 0;
  next = fp.get();
  seq[j] = ' ';
  ++j;
  while((fp.peek()) != EOF){

    if(next == '\n'){
      chk = 1;
    }
    if(next != '\n' && chk == 0){
      name[i] = next;
      ++i;
    }
    if(next != '\n' && chk == 1){
      seq[j] = next;
      ++j;
    }
    next = fp.get();
  }

  seq[seqlength + 1] = '\0';
  fp.close();
}

Direct:: Direct(string &stag, int slen, char *sarr)
// Input parameter stag is a reference to a string object as the name of a sequence.
// Input parameter slen is the length of the sequence.
// Input parameter sarr is a pointer to a char array with the sequence
// saved at index 1 through index slen.

// If slen is negative or sarr is NULL, prints out an error message with fatal. 
// A copy of the string object is saved in the member name.
// The value slen is saved in the member length.
// Memory is allocated to the member seq and a copy of the sequence is saved in seq.
{
  name = stag;
  seq = new char[slen + 2];
  seq = sarr;
}

Direct:: Direct(const Direct &obj)
// A deep copy is made.
{ 
  name = obj.name;
  length = obj.length;
  seq = new char[length + 2];
  seq = obj.seq;
}

Direct:: Direct()
// Sets seq to NULL.
{
  delete[] seq;
  seq = NULL;
}

Direct:: ~Direct()
// Frees heap memory if seq is not NULL.
{
  if(seq != NULL){
    delete[] seq;
    seq = NULL;
  }
}

string Direct:: getName() const
// Returns the member name.
{
  return name; //Do I really need to have this.name or just name?
}

int Direct:: getLength() const
// Returns the member length.
{
  return length;
}

char* Direct:: getSeq() const
// Returns the member seq.
{
  return seq;
}

class Matrix   // a class for an object holding the three matrices
 { private:
    Direct &origin;  // reference to a Direct object called an original sequence
    Direct &derived; // reference to a Direct object called a derived sequence
    int  rowind; // max row index
    int  colind; // max column index
    int  **Dmat; // 2-dimensional array D[rowind + 1][colind + 1]
    int  **Imat; // 2-dimensional array I[rowind + 1][colind + 1]
    int  **Smat; // 2-dimensional array S[rowind + 1][colind + 1]

    int** get2dspace(int rowind, int colind); // allocates a 2D array on heap.
    void free2dspace(int rowind, int **arr);  // frees a 2D array on heap.

    int max(int one, int two);
    int compare(char a, char b, struct scoretp &param);
    int trimax(int one, int two, int three);

   public:
    Matrix(Direct &seqone, Direct &seqtwo, struct scoretp &param); // normal constructor
    ~Matrix();                  // destructor
    Direct& getOrigin() const;  // returns origin.
    Direct& getDerived() const; // returns derived.
    int getRowInd() const;      // returns rowind.
    int getColInd() const;      // returns rowind.
    int **getMat(char kind) const;    // returns a specified matrix.
    string toString(char kind) const; // generates a string form of a specified matrix.
 };

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
  cout<< rowind << "    " << colind <<endl;


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
    max(max(one, two) , three);
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
  int** arr = new int*[rowind + 2];
  for(int i = 0; i < (rowind+2); ++i){
    arr[i] = new int[colind + 2];
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
    
    if(kind == 'S' || kind == 'D' || kind == 'I'){
        for(i = 0; i<=origin.getLength(); ++i){
            for(j=0; j<=derived.getLength(); ++j){
                cout<< getMat(kind)[i][j];
            }
            printf("\n");
        }
    }


}

class Alignment // a class for an object holding an optimal alignment  
 { private:
    Direct &origin;  // reference to a Direct object called an original sequence
    Direct &derived; // reference to a Direct object called a derived sequence
    int  score;      // the score of an optimal alignment of the two sequences
    int  editnum;    // the number of substitutions, insertion and deletion gaps
    int  subinsertlen; // the total length of substitutions and insertion gaps
    int  alignlen;   // the length of the alignment
    char *top;   // a char array of length alignlen for the top row of the alignment
    char *mid;   // a char array of length alignlen for the mid row of the alignment
    char *bot;   // a char array of length alignlen for the bottom row of the alignment

   public:
    Alignment(Matrix &matobj, struct scoretp &param); // normal constructor
    ~Alignment();                                     // destructor
    Direct &getOrigin() const;  // returns origin.
    Direct &getDerived() const; // returns derived.
    int getScore() const; // returns score.
    int getEditNum() const; // returns editnum.
    int getSubInsertLen() const; // returns insertlen.
    int getAlignLen() const; // returns alignlen.
    char* getRow(char kind) const; // returns top, mid, or bot row.
    string toString() const; // generates a string form of the alignment.
 };

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
  int** mat = S;
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
top[k+1] = '\0';
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
  int row = (getAlignLen()/70)+1;

  for(i=0; i<= row; ++i){
        for(j=0; j<=70; ++j){
            if(mid[j] == '\0'){
              topln << top[j];
              fubar << mid[j];
              botln << bot[j];
            }
        }
        fin << topln.str() << '\n';
        fin << fubar.str() << '\n';
        fin << botln.str() << '\n';
    }

  return fin.str();
}

class Encoded        // a class for an object holding compression information.
 { private:
    Direct &origin;  // reference to a Direct object called an original sequence
    char *subinsertion; // an array holding a concatenated sequence of parts in subs and insertions
    int  subinsertlen;  // the length of the concatenated sequence
struct Edit *operation; // an array holding a number of edit operations
    int  editnum;       // the number of edit operations
    string dname;	// the name of the derived sequence
    int   dlength;      // the length of the drived sequence

   public:
    Encoded(Alignment &obj); // a normal constructor
    ~Encoded();              // destructor
    int getEditNum() const;  // returns editnum.
    struct Edit* getOperation() const; // returns operation.
    int getSubInsertLen() const; // returns subinsertlen.
    char* getSubInsertion() const; // return subinsertion.
    int getDLength() const; // return dlength.
    string getDName() const; // return dname.
    Direct& getOrigin() const; // return origin.
    string toString() const; // generates a string form of its contents.
    char* getDSeq() const;   // deirves a sequence and turns it in a char array.
 };

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
  dname = obj.getDerived().getName();
  dlength = obj.getDerived().getLength();
  operation = new struct Edit[editnum+1];
  subinsertion = new char[subinsertlen+1];

  for (int ind = 0; ind < obj.getAlignLen();){
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
  char* A = new char[origin.getLength()];
  char* B = new char[dlength];
  int snind = 0;
  int opind = 0;
  int oppos = 0;
  short indel = 0;
  int k = 0;

  for( ; opind < editnum; opind++){
    oppos = operation[opind].position;
    indel = operation[opind].indel;
    if(i < oppos){
      for(k = 0; k <= oppos-1-i; k++){ // make sure that it is <= or <
        B[j+k] = A[i+k];
      }
      i = oppos;
      j += oppos-i; // need parenthesis?
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
    if(i <= origin.getLength()){
      for(k = 0; k <= origin.getLength()-i; k++){ // make sure that it is <= or <
        B[j+k] = A[i+k];
      }
    }
    B[dlength+1] = '\0';
  }


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

}

int main(int argc, char *argv[])
{
   if ( argc != 6 )
   { cerr << "Usage: " << argv[0] << " Seq1 Seq2 mismatch gap_open gap_extend" << endl << endl;
     cerr << "Seq1        file of one sequence in FASTA format" << endl;
     cerr << "Seq2        file of one sequence in FASTA format" << endl;
     cerr << "mismatch    a negative integer" << endl;
     cerr << "gap_open    gap open penalty, a non-negative integer" << endl;
     cerr << "gap_extend  gap extension penalty, a positive integer" << endl;
     exit(1);
   }

  // Declares and initializes a struct scoretp variable named pararec.
    struct scoretp pararec;
    pararec.match = 10;
    sscanf(argv[3], "%d", &pararec.mismat);
    sscanf(argv[4], "%d", &pararec.gopen);
    sscanf(argv[5], "%d", &pararec.gext);  
  // Declares and constructs a Direct object named seqone for the sequence in file argv[1].
    Direct seqone(argv[1]);
  // Declares and constructs a Direct object named seqtwo for the sequence in file argv[2].
    Direct seqtwo(argv[2]);
  // Declares and constructs a Matrix object named matobj for the matrices.
    Matrix matobj(seqone, seqtwo, pararec);
  // Declares and constructs an Alignment object named alignobj for the alignment.
    Alignment alignobj(matobj, pararec);
    //cout<< alignobj.getRow('T') <<endl;
    //cout<< alignobj.getRow('M') <<endl;
    //cout<< alignobj.getRow('B') <<endl;
    //cout<< matobj.getMat('S')[0][0]<<endl;
    //matobj.toString('S');
  // Declares and constructs an Encoded object named encodobj for the differences.
    Encoded encodobj(alignobj);
  // Prints out the scoring parameters.
    cout<< "Match Score  Mismatch Score  Gap-Open Penalty Gap-Extension Penalty" <<endl;
    cout<< " "<< pararec.match <<"  "<< pararec.mismat <<"  "<< pararec.gopen <<"  "<< pararec.gext <<endl;
  // Sends the alignment along with other information to cout by calling toString() from alignobj.
   cout << alignobj.toString() << endl;
  // Sends the contents in the object encodobj to cout by calling toString() from it.
    cout<< encodobj.toString() << endl;
  // Gets the sequence in seqtwo by calling getSeq().
    char* st = seqtwo.getSeq();
  // Derives the sequence in dseq by calling getDSeq() from encodobj.
    char* dse = encodobj.getDSeq();
  // Runs strcmp() on the two sequences and prints out one of the two
  // messages based on the comparison:
  //  The original and derived sequences are identical.
  //  The original and derived sequences are different.
    if(strcmp(st,dse) == 0){
      cout<<"The original and derived sequences are identical."<<endl;
    }else{
      cout<<"The original and derived sequences are different."<<endl;
    }
  // Prints out each of the two sequences.
    cout<<"Original: "<< st <<endl;
    cout<<"Derived: "<< dse <<endl;
  // sends each matrix to cout by calling the toString() of Matrix from matobj.
    cout<<matobj.toString('S')<<endl;
    cout<<matobj.toString('D')<<endl;
    cout<<matobj.toString('I')<<endl;

  return 0;
}

// Prints out an error message and stops.
void fatal(const char *msg)
{
     cerr << "Error message: " << string(msg) << endl;
     exit(1); // stops unexpectedly.
}

// Prints out an error message and stops.
void fatal(const char *msg, const char *name) 
{
     cerr << "Error message: " << string(msg) << string(name) << endl;
     exit(1); // stops unexpectedly.
}

// Opens an input file stream and checks for success.
void ckopeninf(ifstream &infile, const char *fname)
{
	infile.open(fname);
	if ( infile.fail() )
	  fatal("ckopeninf: cannot open input file ", fname);
}

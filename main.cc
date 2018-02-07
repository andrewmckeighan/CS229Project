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
#include "template.stack.min.hh"
#include "Encoded.hh"
#include "prototype.hh"

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
    pararec.mismat = -20;
    pararec.gopen = 40;
    pararec.gext = 2;
  // Declares and constructs a Direct object named seqone for the sequence in file argv[1].
    Direct origin(argv[1]);
  // Declares and constructs a Direct object named seqtwo for the sequence in file argv[2].
    Stack<Encoded>* enc = new Stack<Encoded>;
    Stack<Compressed>* com = new Stack<Compressed>;

    int track = 2;
    while(track <= argc-1){
      Direct derived(argv[track]);
      Matrix matobj(origin, derived, pararec);
      Alignment alignobj(matobj, pararec);
      //Encoded encodobj(alignobj);
      Encoded* ep = new Encoded(alignobj);
      //Compressed compobj(alignobj);
      Compressed* cp = new Compressed(alignobj);
      enc->push(ep);
      com->push(cp);
      //delete ep;
      //delete cp;
      track++;
    }
    Encoded* emin = findMin(enc);
    Encoded* emax = findMax(enc);
    Compressed* cmin = findMin(com);
    Compressed* cmax = findMax(com);

    char* emaa = emax->getDSeq();
    string eman = emax->getDName();

    char* ema = emin->getDSeq();
    string emn = emin->getDName();

    char* cx = cmax->getDSeq();
    string cxn = cmax->getDName();
    char* cm = cmin->getDSeq();
    string cmn = cmin->getDName();
    


    cout<< "\nA stack of Encoded pointers"<<endl;
    cout<< "Min's number of differences: "<< emin->getNumDiff() <<endl;//not done!
    cout<< "Min's name : "<< emn <<endl;
    cout<< "Min's derived sequence: "<< ema <<"\n"<<endl;

//    james = std::str(emaa);
    cout<< "Max's number of differences: "<< emax->getNumDiff() <<endl;//not done!
    cout<< "Max's name : "<< eman <<endl;
    cout<< "Max's derived sequence: "<< emaa <<"\n"<<endl;

  //  james = std::str(cm);
    cout<< "A stack of Compressed pointers"<<endl;
    cout<< "Min's number of differences: "<< cmin->getNumDiff() <<endl; //not done!
    cout<< "Min's name : "<< cmn <<endl;
    cout<< "Min's derived sequence: "<< cm <<"\n"<<endl;

//    std::string str(cx);
    cout<< "Max's number of differences: "<< cmax->getNumDiff() <<endl;//not done!
    cout<< "Max's name : "<< cxn <<endl;
    cout<< "Max's derived sequence: "<< cx <<"\n"<<endl;


  // Sends the alignment along with other information to cout by calling toString() from alignobj.
  //  cout << alignobj.toString() << endl;
  // Sends the contents in the object encodobj to cout by calling toString() from it.
  //  cout<< encodobj.toString() << endl;
  // Gets the sequence in seqtwo by calling getSeq().
  //  char* st = seqtwo.getSeq();
  // Derives the sequence in dseq by calling getDSeq() from encodobj.
  //  char* dse = encodobj.getDSeq();
  // Free everthang yo!

  return 0;
}

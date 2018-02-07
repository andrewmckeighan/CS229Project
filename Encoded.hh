//@author Andrew McKeighan
#ifndef Encoded_G_H  // A header guard is introduced to avoid including
#define Encoded_G_H  // this header file twice in any source code file.

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
    bool operator<=(Encoded &rightobj) const; //to compare two Encoded objects
	int getNumDiff() const; 
 };

class Compressed : public Encoded
 { 
   public:
    Compressed(Alignment &obj);
    bool operator<=(Encoded &rightobj) const;
};

#endif
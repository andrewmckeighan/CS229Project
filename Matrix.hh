//@author Andrew McKeighan
#ifndef Matrix_G_H  // A header guard is introduced to avoid including
#define Matrix_G_H  // this header file twice in any source code file.


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
#endif
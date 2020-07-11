#ifndef ADI2D_HPP
#define ADI2D_HPP

#include <suitesparse/umfpack.h>
#include <iostream>

// Held const -- means that changes in fcn won't matter
static const uint16_t nrows = 21;
static const uint16_t ncols = 21;

// Temperature: the array that gets modified over time
static double T[nrows][ncols];

class Adi{

public:
    float t = 0; // Current time
    float t_max = 10.; // Time at end of run
    float D = 1.; // Diffusivity
    float dx = 1.; // Spatial step
    float dt = 1.; // Time step
    double K; // Lumped coefficient: diffusivity, space, time

    // Nonzero rows, top to bottom, for each column in the
    // sparse matrix A
    int Ai_rows[3*ncols-2];
    int Ai_cols[3*nrows-2];

    // Denotes which entries within Ai belong to which
    // column, in order. Inclusive of first index,
    // exclusive of second, incrememnts by 1 for each row.
    int Ap_rows[ncols+1];
    int Ap_cols[nrows+1];

    // Numerical values in sparse matrix.
    // These are ordered first by column and then by row.
    double Ax_rows[3*ncols-2];
    double Ax_cols[3*nrows-2];

    // UMFPACK pointers
    void *Symbolic;
    void *Numeric;

    // Constructor
    Adi();

    // Functions
    void initialize();
    void update();
    void run();
    void finalize();

private:
    void _set_A_matrix_row( uint16_t _row_index );
    void _set_A_matrix_column( uint16_t _column_index );
    void _update_rows();
    void _update_columns();

};

#endif

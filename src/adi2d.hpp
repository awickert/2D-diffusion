#ifndef ADI2D_HPP
#define ADI2D_HPP

#include <umfpack>
#include <iostream>

class adi{

public:
    uint8_t nrows;
    uint8_t ncols;
    float K; // Lumped coefficient: diffusivity, space, time
    float t = 0; // Current time

    // Nonzero rows, top to bottom, for each column in the
    // sparse matrix A
    uint16_t *Ai_rows[3*ncols-2]; 
    uint16_t *Ai_cols[3*nrows-2];
 
    // Denotes which entries within Ai belong to which
    // column, in order. Inclusive of first index,
    // exclusive of second, incrememnts by 1 for each row.
    uint16_t *Ap_rows[ncols+1];
    uint16_t *Ap_cols[nrows+1];

    // Numerical values in sparse matrix.
    // These are ordered first by column and then by row.
    float *Ax_rows[3*ncols-2];
    float *Ax_cols[3*nrows-2];
    
    // UMFPACK pointers
    void *Symbolic;
    void *Numeric;

    void initialize(uint8_t nrows = 21, uint8_t ncols = 21, float D = 1.,
                    float dx = 1., float dt=1., float t_max=10.);
    void update();
    void run();
    void finalize();

private:
    void set_A_matrix_row( uint8_t _row_index );
    void set_A_matrix_column( uint8_t _column_index );
    void _update_rows();
    void _update_columns();

#endif

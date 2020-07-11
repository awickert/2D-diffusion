#include "adi2d.hpp"

void initialize(uint8_t nrows, uint8_t ncols, float D, float dx, float dt, float t_max){
    // Coefficient (constant)
    K = D * dt / (dx*dx)
    // T array
    float T[nrows][ncols];
    for (int i = 0; i<ncols; i++){
        for (int j = 0; j<nrows; j++){
            T[i,j] = 0.;
        }
    }
    // Initial high T in the center
    T[i/2,j/2] = 1.

    // Array A for each row within array T
    set_A_matrix_row()
    // and for each column
    //set_A_matrix_column()
}

void set_A_matrix_row( uint8_t _row_index ){
    // Rows within each column, in order
    // Two in first column
    Ai_rows[0] = 0;
    Ai_rows[1] = 1;
    // Three for all intermediary ones
    for (int i = 1; i<(ncols-1); i++){
        Ai_rows[3*i-1] = i-1; // upper
        Ai_rows[3*i] = i; // main diagonal
        Ai_rows[3*i+1] = i+1; // lower
    }
    // Two at the end
    Ai_rows[3*ncols-4] = ncols-2;
    Ai_rows[3*ncols-3] = ncols-1;

    // Pointers to rows in Ai within each column, in order
    // Two rows in the first and last columns, three in the middle
    // Indices are (inclusive, exclusive)
    Ap_rows[0] = 0;
    Ap_rows[1] = 2;
    for (int i = 2; i<(ncols-1); i++){
        Ap_rows[i] = Ap_rows[i-1] + 3;
    Ap_rows[ncols] = Ap_rows[ncols-1] + 2;

    // And finally, we set the values within each of these rows.
    // This solution is for uniform diffusivity, hence this need be done
    // only once. For TWSM, these values will need to be updated every time
    // (Hence my passing _row_index but not using it)
    // step, for every row and column.
    // Two elements in the first column
    Ax_rows[0] = 2*K+1;
    Ai_rows[1] = -1*K;
    // Three for all intermediary ones
    for (int i = 1; i<(ncols-1); i++){
        Ai_rows[3*i-1] = -1*K; // upper
        Ai_rows[3*i] = 2*K+1; // main diagonal
        Ai_rows[3*i+1] = -1*K; // lower
    }
    // Two at the end
    Ai_rows[3*ncols-4] = -1*K;
    Ai_rows[3*ncols-3] = 2*K+1;
}

void _update_rows(){
    for i in range(nrows){
        float* _T = T[i]
        umfpack_di_symbolic( ncols, ncols, Ap_rows, Ai_rows, Ax_rows,
                             &Symbolic,
                             null, null );
        umfpack_di_numeric( Ap_rows, Ai_rows, Ax_rows,
                            Symbolic, &Numeric,
                            null, null );
        umfpack_di_free_symbolic( &Symbolic );
        umfpack_di_solve( UMFPACK_A, Ap_rows, Ai_rows, Ax_rows,
                            T[i], T[i],
                            Numeric,
                            null, null )
        umfpack_di_free_numeric (&Numeric);
    }
}

void set_A_matrix_column( uint8_t _column_index ){

void _update_columns(){

}

void update(){
    update_rows();
    update_columns();
}

void run(){
    while( t < t_max ){
        update();
    }
}

void finalize(){
    cout << T[i/2,j/2] << ".\n";
}

int main(){
  cout << "Test!" << ".\n";
}

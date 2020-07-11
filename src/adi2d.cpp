#include "adi2d.hpp"

Adi::Adi(){}

void Adi::initialize(uint16_t nrows, uint16_t ncols, float D, float dx, float dt, float t_max){
    // Number of rows and columns
    nrows = nrows;
    ncols = ncols;
    // Coefficient (constant); treated here as if dx = dy
    // otherwise would have two "K" variables
    K = D * dt / (dx*dx);
    // Time step
    dt = dt;
    // Maximum time at end of run
    t_max = t_max;
    // T array
    for (uint16_t i = 0; i<ncols; i++){
        for (uint16_t j = 0; j<nrows; j++){
            T[i][j] = 0.;
        }
    }
    // Initial high T in the center
    T[nrows/2][ncols/2] = 1.;

    // Array A for each row within array T
    set_A_matrix_row( 1 ); // Dummy index for now
    // and for each column
    //set_A_matrix_column()
}

void Adi::set_A_matrix_row( uint16_t _row_index ){
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
    for (int i = 2; i<(ncols); i++){
        Ap_rows[i] = Ap_rows[i-1] + 3;
    }
    Ap_rows[ncols] = Ap_rows[ncols-1] + 2;

    // And finally, we set the values within each of these rows.
    // This solution is for uniform diffusivity, hence this need be done
    // only once. For TWSM, these values will need to be updated every time
    // (Hence my passing _row_index but not using it)
    // step, for every row and column.
    // Two elements in the first column
    Ax_rows[0] = 2*K+1;
    Ax_rows[1] = -1*K;
    // Three for all intermediary ones
    for (uint32_t i = 1; i<(ncols-1); i++){
        Ax_rows[3*i-1] = -1*K; // upper
        Ax_rows[3*i] = 2*K+1; // main diagonal
        Ax_rows[3*i+1] = -1*K; // lower
    }
    // Two at the end
    Ax_rows[3*ncols-4] = -1*K;
    Ax_rows[3*ncols-3] = 2*K+1;
}

void Adi::_update_rows(){
    for (uint16_t i = 0; i<nrows; i++){
        double* _T = T[i];
        umfpack_di_symbolic( ncols, ncols, Ap_rows, Ai_rows, Ax_rows,
                             &Symbolic,
                             NULL, NULL );
        umfpack_di_numeric( Ap_rows, Ai_rows, Ax_rows,
                            Symbolic, &Numeric,
                            NULL, NULL );
        umfpack_di_free_symbolic( &Symbolic );
        umfpack_di_solve( UMFPACK_A, Ap_rows, Ai_rows, Ax_rows,
                            _T, T[i],
                            Numeric,
                            NULL, NULL );
        //for (i = 0 ; i < ncols ; i++) printf ("x [%d] = %g\n", i, _T [i]) ;
        for (int j=0; j<(ncols); j++){
            //std::cout << T[i][j] << " ";
            std::cout << _T[j] << " ";
            //std::cout << Ap_rows[j] << " ";
        }
        std::cout << "\n";
        //std::cout << "Test!" << ".\n";
        //umfpack_di_free_numeric (&Numeric);
    }
}

void Adi::set_A_matrix_column( uint16_t _column_index ){

}

void Adi::_update_columns(){

}

void Adi::update(){
    _update_rows();
    _update_columns();
}

void Adi::run(){
    while( t < t_max ){
        update();
        t += dt;
    }
}

void Adi::finalize(){
    //std::cout << T[nrows/2][ncols/2] << "\n";
    //std::cout << K << "\n";
}

int main(){
    std::cout << "Test!" << "\n";
    Adi adi;
    adi.initialize();
    adi.update();
    adi.finalize();
    //std::cout << T[nrows/2][ncols/2] << ".\n";
}

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>

#include <gsl/gsl_math.h>
#include <gsl/gsl_eigen.h>

#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>

#include <gsl/gsl_linalg.h>
#include <gsl/gsl_blas.h>



int main(void)
{
    int i, j;

    //<summary> This part is getting the matrix x_train real part from matlab
    // the matrix from matlab we divided it into real part and image part
    // because of  gsl loading data's way,  we need to transpose the matrix before we load in
    // for example: the data at Matlab is 640*1062 we need to transpose the matrix to 1062*640
    // and we claim that a register called 640*1062 , the data type is default double  
    gsl_matrix* x_train_real = gsl_matrix_alloc(640, 1062);

    {
        FILE* file_x_train_real = fopen("x_train_real_tp.bin", "rb");
        gsl_matrix_fread(file_x_train_real, x_train_real);
        fclose(file_x_train_real);
    }


    /*  
    * this is the way we test about load in matrix,the matrix is 640*1062
    * 
    for (i = 0; i < 640; i++)
    {
        for (j = 0; j <1062; j++)
        {
            double real_number = gsl_matrix_get(x_train_real, i, j);
            printf("this is the element %d ,%d = %g\n",i,j, real_number);
        }

    }
    */


    ///<summary> This part is getting the matrix x_train image part from matlab

    gsl_matrix* x_train_imag = gsl_matrix_alloc(640, 1062);

    {
        FILE* file_x_train_imag = fopen("x_train_imag_tp.bin", "rb");
        gsl_matrix_fread(file_x_train_imag, x_train_imag);
        fclose(file_x_train_imag);
    }

    /* this part is testing for the matrix x_train image part ,we load from Matlab 
    * 
    i = 315;
    j = 959;

    printf("this is the element %d ,%d = %g\n", i, j, gsl_matrix_get(x_train_imag, i, j));

    */


    /// <summary> x_train_complex
    // this part is combine the matrix x_train real part and x_train image part together
    // and we store the complex data into the complex matrix 
    gsl_matrix_complex* x_train_complex = gsl_matrix_complex_alloc(640, 1062);
    for (i = 0; i < 640; i++)
    {
        for (j = 0; j < 1062; j++)
        {
            double x_train_real_part = gsl_matrix_get(x_train_real, i, j);
            double x_train_imag_part = gsl_matrix_get(x_train_imag, i, j);

            gsl_matrix_complex_set(x_train_complex, i, j, gsl_complex_rect(x_train_real_part, x_train_imag_part));
        }


    }

    /* this is the testing part of x_train_complex, we want to know whethere it is combine correctly
    
    i = 122;
    j = 827;

    gsl_complex x_train_complex_i = gsl_matrix_complex_get(x_train_complex, i,j);
    printf("this is the element of x_train_complex %d ,%d = %g + %g i\n", i, j, GSL_REAL(x_train_complex_i), GSL_IMAG(x_train_complex_i));
    */



    /// <summary> this part is about loading vector sv_real ,the real part of vector "sv".
    // we load the data from matlab , but we dont need to transpose it, because it is just a vector
    // we claim that the vector is 640 
    gsl_vector* sv_real = gsl_vector_alloc(640);

    {
        FILE* file_sv_real = fopen("sv_real.bin", "rb");
        gsl_vector_fread(file_sv_real, sv_real);
        fclose(file_sv_real);
    }
    // this is the testing part of sv_real vector , whethere it is load in correctly
    /*
    for (i = 0; i < 640; i++)
    {
        printf("%g\n", gsl_vector_get(sv_real, i));
    }
    */

    /// <summary> this part is about loading vector sv_imag ,the image part of vector "sv".
    // we load the data from matlab , but we dont need to transpose it, because it is just a vector
    // we claim that the vector is 640 
    gsl_vector* sv_imag = gsl_vector_alloc(640);

    {
        FILE* file_sv_imag = fopen("sv_imag.bin", "rb");
        gsl_vector_fread(file_sv_imag, sv_imag);
        fclose(file_sv_imag);
    }

    /* this is the testing part of sv_image vector , whethere it is load in correctly
    for (i = 0; i < 640; i++)
    {
        printf("%g\n", gsl_vector_get(sv_imag, i));
    }
    */


    ////<summary> we combine the vector about real part and image part , the real and image is setting into the sv_complex.
    gsl_vector_complex* sv_complex = gsl_vector_complex_alloc(640);
    for (j = 0; j < 640; j++)
    {
        double sv_real_part = gsl_vector_get(sv_real, j);
        double sv_imag_part = gsl_vector_get(sv_imag, j);

        gsl_vector_complex_set(sv_complex, j, gsl_complex_rect(sv_real_part, sv_imag_part));
    }

    /* this is the part we are testing the sv_complex , whethere it is combine correctly
     i = 123;
        gsl_complex sv_complex_i = gsl_vector_complex_get(sv_complex, i);

        printf("the complex vector sv i(%d) is  =  %g + %gi\n", i, GSL_REAL(sv_complex_i), GSL_IMAG(sv_complex_i));

    */

    /// this is the all text , we test about , whethere the sv_complex and the x_train_complex is correct 
    /*
    i = 123;
    j = 456;
    gsl_complex sv_complex_i = gsl_vector_complex_get(sv_complex, i);
    printf("the complex vector sv i(%d) is  =  %g + %gi\n", i, GSL_REAL(sv_complex_i), GSL_IMAG(sv_complex_i));

    gsl_complex x_train_complex_i = gsl_matrix_complex_get(x_train_complex, i, j);
    printf("this is the element of x_train_complex %d ,%d = %g + %g i\n", i, j, GSL_REAL(x_train_complex_i), GSL_IMAG(x_train_complex_i));

    */

    gsl_vector_complex* x_train_get_col = gsl_vector_complex_alloc(640);
    gsl_matrix_complex* x_train_get_col_matrix = gsl_matrix_complex_alloc(640, 1);
    gsl_matrix_complex* x_train_get_col_matrix_ans = gsl_matrix_complex_alloc(640, 640);
    gsl_matrix_complex* Ru = gsl_matrix_complex_calloc(640, 640);
    gsl_matrix_complex* Ru_inverse = gsl_matrix_complex_calloc(640, 640);


    int idex;

    for (idex = 0; idex < 1062; idex++)
    {
        // we get the matrix "x_train_complex" 's one column, and we store tthat column in the vector "x_train_get_col"
        gsl_matrix_complex_get_col(x_train_get_col, x_train_complex, idex);

        // we store the  vector "x_train_get_col" into a "640*1" matrix "x_train_get_col_matrix"
        gsl_matrix_complex_set_col(x_train_get_col_matrix, 0, x_train_get_col);

        //we use this function to do the matrix * matrix. which means we want to get the "A * A ' " , Matrix multiple Conjugate transpose (Matrix).
        // and we store the A * A' into "x_train_get_col_matrix_ans"  
        gsl_blas_zgemm(CblasNoTrans, CblasConjTrans, gsl_complex_polar(1.0, 0.0), x_train_get_col_matrix, x_train_get_col_matrix, gsl_complex_polar(0.0, 0.0), x_train_get_col_matrix_ans);
        
        //we use this function to do the equation : Ru = Ru + A* A' ;
        gsl_matrix_complex_add(Ru, x_train_get_col_matrix_ans);
        printf("success %d \n", idex);
    }

    /*  this is the part we test about the   "matrix * matrix' + Ru = Ru "
    i = 303;
    j = 3;
    gsl_complex x_train_complex_ans_i = gsl_matrix_complex_get(Ru, i, j);
    //printf("this is the element of matrix*matrix' col = 0 is  %d ,%d = %g + %g i\n", i, j, GSL_REAL(x_train_complex_ans_i), GSL_IMAG(x_train_complex_ans_i));
    printf("this is the element of matrix*matrix' + Ru = Ru is  %d ,%d = %g + %g i\n", i, j, GSL_REAL(x_train_complex_ans_i), GSL_IMAG(x_train_complex_ans_i));
    */

    // this part is to spawn the Ru for testing the matrix form 
    {
        FILE* file_Ru_formatlab = fopen("Ru_for_matlab.bin", "wb");
        gsl_matrix_complex_fwrite(file_Ru_formatlab, Ru);
        fclose(file_Ru_formatlab);
    }


    gsl_permutation* p = gsl_permutation_alloc(640);
    int s;
    gsl_linalg_complex_LU_decomp(Ru, p, &s);

    /* this is the part we are testing the RU decompose
    i = 0;
    j = 0;
    gsl_complex x_train_complex_ans_iii = gsl_matrix_complex_get(Ru, i, j);
    printf("this is the element of RU decompose is  %d ,%d = %g + %g i\n", i, j, GSL_REAL(x_train_complex_ans_iii), GSL_IMAG(x_train_complex_ans_iii));
    */

    gsl_linalg_complex_LU_invert(Ru, p, Ru_inverse);

    /*
    //test for all  Ru_inverse
    i = 0;
    j = 0;
    gsl_complex x_train_complex_ans_i = gsl_matrix_complex_get(Ru_inverse, i, j);
    printf("this is the element of Ru_inverse = Ru is  %d ,%d = %g + %g i\n", i, j, GSL_REAL(x_train_complex_ans_i), GSL_IMAG(x_train_complex_ans_i));
    */

    // this part is to spawn the Ru inverse for testing the matrix form 
    {
        FILE* file_Ru_inverse_formatlab = fopen("Ru_inverse_for_matlab.bin", "wb");
        gsl_matrix_complex_fwrite(file_Ru_inverse_formatlab, Ru_inverse);
        fclose(file_Ru_inverse_formatlab);
    }


    gsl_vector_complex* Ru_inv_mul_sv = gsl_vector_complex_alloc(640);

    // we use this function to do the matrix * vector ,which the equation is " inv(Ru) * sv_complex ", and we store the answer vector into the Ru_inv_mul_sv
    gsl_blas_zgemv(CblasNoTrans, GSL_COMPLEX_ONE, Ru_inverse, sv_complex, GSL_COMPLEX_ZERO, Ru_inv_mul_sv);


    // this part is we are test the Ru_inv_mul_sv , which it is the weight vector , the answer we want .
    i = 127;
    gsl_complex w_temp = gsl_vector_complex_get(Ru_inv_mul_sv, i);
    printf("the complex vector w_temp i(%d) is  =  %g + %gi\n", i, GSL_REAL(w_temp), GSL_IMAG(w_temp));


    // this part is we store the vector "Ru_inv_mul_sv" into the matrix (640*1) "Ru_inv_mul_sv_matrix"
    gsl_matrix_complex* Ru_inv_mul_sv_matrix = gsl_matrix_complex_alloc(640, 1);
    gsl_matrix_complex_set_col(Ru_inv_mul_sv_matrix, 0, Ru_inv_mul_sv);


    // this part is we write the matrix Ru_inv_mul_sv_matrix (answer) into the  file "w_temp_matrix.bin" , which it is writing in binary 
    // must remember these data dont need to transpose angain, we can use it directly in matlab 
    // and the data in matalb if we want to use it, we need to know 
    // the data is sorting in 
    // real part of first element
    // image part of first element
    // real part of second element
    // image part of second element
    //.....

    {
        FILE* file_w_temp = fopen("w_temp_matrix.bin", "wb");
        gsl_matrix_complex_fwrite(file_w_temp, Ru_inv_mul_sv_matrix);
        fclose(file_w_temp);
    }


    gsl_matrix_complex_free(Ru_inv_mul_sv_matrix);
    gsl_vector_complex_free(Ru_inv_mul_sv);
    gsl_vector_complex_free(x_train_get_col);
    gsl_matrix_complex_free(x_train_get_col_matrix);
    gsl_matrix_complex_free(x_train_get_col_matrix_ans);
    gsl_matrix_complex_free(Ru);
    gsl_matrix_complex_free(Ru_inverse);

    gsl_vector_free(sv_real);
    gsl_vector_free(sv_imag);
    gsl_vector_complex_free(sv_complex);

    gsl_matrix_free(x_train_real);
    gsl_matrix_free(x_train_imag);
    gsl_matrix_complex_free(x_train_complex);

    return 0;
}
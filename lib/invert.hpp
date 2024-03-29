#ifndef _INVERT_HPP
#define _INVERT_HPP

// #include <cstdio>
// #include <vector>

// // https://github.com/atchekho/harmpi/blob/master/lu.c

// /*************************************************************************/
// /*************************************************************************
//    LU_decompose():
//        Performs a LU decomposition of the matrix A using Crout's method      
//        with partial implicit pivoting.  The exact LU decomposition of the    
//        matrix can be reconstructed from the resultant row-permuted form via  
//        the integer array permute[]                                            
                                                                             
//        The algorithm closely follows ludcmp.c of "Numerical Recipes  
//        in C" by Press et al. 1992.                                           
                                                                             
//        This will be used to solve the linear system  A.x = B                 
                                                                             
//        Returns (1) if a singular matrix is found,  (0) otherwise.            
// *************************************************************************/
// template <typename M,
// 	  typename T=typename M::type>
// int LU_decompose(M& A,std::vector<int>& permute)
// {
//   const int nR=A.nR;

//   std::vector<T> row_norm(nR);
//   T absmax,maxtemp,mintemp;
  
//   int i, j, k, max_row;


//    max_row = 0;

//   /* Find the maximum elements per row so that we can pretend later
//      we have unit-normalized each equation: */

//   for( i = 0; i < nR; i++ ) { 
//     absmax = 0.;
    
//     for( j = 0; j < nR ; j++ ) { 
      
//       maxtemp = abs( A(i,j) ); 

//       if( maxtemp > absmax ) { 
// 	absmax = maxtemp; 
//       }
//     }

//     /* Make sure that there is at least one non-zero element in this row: */
//     if( absmax == 0. ) { 
//      fprintf(stderr, "LU_decompose(): row-wise singular matrix!\n");
//       return(1);
//     }

//     row_norm[i] = 1. / absmax ;   /* Set the row's normalization factor. */
//   }


//   /* The following the calculates the matrix composed of the sum 
//      of the lower (L) tridagonal matrix and the upper (U) tridagonal
//      matrix that, when multiplied, form the original maxtrix.  
//      This is what we call the LU decomposition of the maxtrix. 
//      It does this by a recursive procedure, starting from the 
//      upper-left, proceding down the column, and then to the next
//      column to the right.  The decomposition can be done in place 
//      since element {i,j} require only those elements with {<=i,<=j} 
//      which have already been computed.
//      See pg. 43-46 of "Num. Rec." for a more thorough description. 
//   */

//   /* For each of the columns, starting from the left ... */
//   for( j = 0; j < nR; j++ ) {

//     /* For each of the rows starting from the top.... */

//     /* Calculate the Upper part of the matrix:  i < j :   */
//     for( i = 0; i < j; i++ ) {
//       for( k = 0; k < i; k++ ) { 
// 	A(i,j) -= A(i,k) * A(k,j);
//       }
//     }

//     absmax = 0.0;

//     /* Calculate the Lower part of the matrix:  i <= j :   */

//     for( i = j; i < nR; i++ ) {

//       for (k = 0; k < j; k++) { 
// 	A(i,j) -= A(i,k) * A(k,j);
//       }

//       /* Find the maximum element in the column given the implicit 
// 	 unit-normalization (represented by row_norm[i]) of each row: 
//       */
//       maxtemp = abs(A(i,j)) * row_norm[i] ;

//       // if( maxtemp >= 1absmax ) {
// 	absmax = maxtemp;
// 	max_row = i;
//       // }

//     }

//     /* Swap the row with the largest element (of column j) with row_j.  absmax
//        This is the partial pivoting procedure that ensures we don't divide
//        by 0 (or a small number) when we solve the linear system.  
//        Also, since the procedure starts from left-right/top-bottom, 
//        the pivot values are chosen from a pool involving all the elements 
//        of column_j  in rows beneath row_j.  This ensures that 
//        a row  is not permuted twice, which would mess things up. 
//     */
//     if( max_row != j ) {

//       /* Don't swap if it will send a 0 to the last diagonal position. 
// 	 Note that the last column cannot pivot with any other row, 
// 	 so this is the last chance to ensure that the last two 
// 	 columns have non-zero diagonal elements.
//        */

//       if( (j == (nR-2)) && (A(j,j+1) == 0.) ) {
// 	max_row = j;
//       }
//       else { 
// 	for( k = 0; k < nR; k++ ) { 

// 	  maxtemp       = A(   j   ,k) ; 
// 	  A(   j   ,k) = A(max_row,k) ;
// 	  A(max_row,k) = maxtemp; 

// 	}

// 	/* Don't forget to swap the normalization factors, too... 
// 	   but we don't need the jth element any longer since we 
// 	   only look at rows beneath j from here on out. 
// 	*/
// 	row_norm[max_row] = row_norm[j] ; 
//       }
//     }

//     /* Set the permutation record s.t. the j^th element equals the 
//        index of the row swapped with the j^th row.  Note that since 
//        this is being done in successive columns, the permutation
//        vector records the successive permutations and therefore
//        index of permute[] also indexes the chronology of the 
//        permutations.  E.g. permute[2] = {2,1} is an identity 
//        permutation, which cannot happen here though. 
//     */

//     permute[j] = max_row;

//     // if( A(j,j) == 0. ) { 
//     //   A(j,j) = absmin;
//     // }


//   /* Normalize the columns of the Lower tridiagonal part by their respective 
//      diagonal element.  This is not done in the Upper part because the 
//      Lower part's diagonal elements were set to 1, which can be done w/o 
//      any loss of generality.
//   */
//     if( j != (nR-1) ) { 
//       maxtemp = 1. / A(j,j)  ;
      
//       for( i = (j+1) ; i < nR; i++ ) {
// 	A(i,j) *= maxtemp;
//       }
//     }

//   }

//   return(0);

//   /* End of LU_decompose() */

// }

// /************************************************************************
//    LU_substitution():
//        Performs the forward (w/ the Lower) and backward (w/ the Upper)   
//        substitutions using the LU-decomposed matrix A[][] of the original
//        matrix A' of the linear equation:  A'.x = B.  Upon entry, A[][]   
//        is the LU matrix, B[] is the source vector, and permute[] is the  
//        array containing order of permutations taken to the rows of the LU
//        matrix.  See LU_decompose() for further details. 		     
    								     
//        Upon exit, B[] contains the solution x[], A[][] is left unchanged.
								     
// ************************************************************************/

// template <typename M,
// 	  typename T=typename M::type>
// void LU_substitution( M& A,std::vector<T>& B,std::vector<int>& permute )
// {
//   int i, j ;
//   int nR = A.nR;
//   T tmpvar,tmpvar2;

  
//   /* Perform the forward substitution using the LU matrix. 
//    */
//   for(i = 0; i < nR; i++) {

//     /* Before doing the substitution, we must first permute the 
//        B vector to match the permutation of the LU matrix. 
//        Since only the rows above the currrent one matter for 
//        this row, we can permute one at a time. 
//     */
//     tmpvar        = B[permute[i]];
//     B[permute[i]] = B[    i     ];
//     for( j = (i-1); j >= 0 ; j-- ) { 
//       tmpvar -=  A(i,j) * B[j];
//     }
//     B[i] = tmpvar; 
//   }
	   

//   /* Perform the backward substitution using the LU matrix. 
//    */
//   for( i = (nR-1); i >= 0; i-- ) { 
//     for( j = (i+1); j < nR ; j++ ) { 
//       B[i] -=  A(i,j) * B[j];
//     }
//     B[i] /= A(i,i) ; 
//   }

//   /* End of LU_substitution() */
// }

// template <typename M,
// 	  typename T=typename M::type>
// int invert_matrix(M Amtmp,M& Aminv)
// {
//   using namespace std;
  
//   const int nR=Amtmp.nR;
//   vector<int> permute(nR);
//   vector<T> dxm(nR);
  
//   // Get the LU matrix:
//   if(LU_decompose(Amtmp,permute)!=0)
//     {
//       fprintf(stderr, "invert_matrix(): singular matrix encountered! \n");
//       return(1);
//     }
  
//   for(int i=0;i<nR;i++)
//     {
//       for(int j=0;j<nR;j++)
// 	dxm[j]=0;
//       dxm[i]=1;
      
//     /* Solve the linear system for the i^th column of the inverse matrix: :  */
//     LU_substitution( Amtmp,  dxm, permute );
    
//     for(int j=0;j<nR;j++)
//       Aminv(j,i)=dxm[j];
//     }
  
//   return(0);
// }

#endif

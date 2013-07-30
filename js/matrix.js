// Returns a matrix that is defined by str. Each row of the matrix is separated
// by `row`, each element of the row is separated by `col`.
function matrix_make_from_string(str,row,col){
	
	// Split the input into rows
	var rows = str.split(row);
	
	// Create a blank matrix of empty rows
	var h = rows.length;
	var M = [];
	var i=0, j=0;
	for(i=0; i<h; i+=1){ M[i] = []; }
	
	// Iterate over each row and then insert each element into that row
	var cols = []; var w = 0;
	for(i=0; i<h; i+=1){
		cols = rows[i].split(col);
		w = cols.length; // this w should NOT change between rows
		for(j=0; j<w; j+=1){
			M[i][j] = parseFloat(cols[j]); // maybe we want a parseInt
		}
	}
	
	// Return the matrix
	return M;
}

// Returns a `w` by `h` matrix of 0's (`h` rows and `w` columns)
function matrix_make_zero(w,h){
	// Create the blank rows
	var R = [];
	var i=0, j=0;
	for(i=0; i<h; i+=1){ R[i] = []; }
	
	// Insert 0's everywhere
	for(i=0; i<h; i+=1){
		for(j=0; j<w; j+=1){
			R[i][j] = 0;
		}
	}
	
	// Return the matrix
	return R;
}

// Returns an `n` by `n` identity matrix (1's on diagonal, 0's elsewhere)
function matrix_make_identity(n){
	var i=0, j=0, I = [];
	for(i=0; i<n; i+=1){
		I[I.length]=[];
		for(j=0; j<n; j+=1){
			// if we're on the diagonal, put a 1; otherwise a 0
			if(i===j){ I[i][j] = 1; }
			else{     I[i][j] = 0; }
		}
	}
	return I;
}

// Returns the trace (scalar) of a matrix. The trace is the sum of the
// elements on the diagonal.
function matrix_trace(M){
	// Make sure the matrix is square (error)
	var h = M.length;
	var w = M[0].length;
	if( h !== w ){ return [[]]; }

	// Sum the diagonal
	var i, sum = 0;
	for( i=0; i<w; i+=1 ){
		sum += M[i][i];
	}
	
	// Return the sum
	return sum;
}

// Returns the psuedo-inverse of matrix M. The psuedo inverse is defined as
// `(AT*A)'*AT` and is also referred to as 'dagger' 'conjugate transpose' 
function matrix_psuedo_inverse(M){
	var MT = matrix_transpose(M);
	var MTM = matrix_multiply(MT,M);
	MTM = matrix_inverse(MTM);
	return matrix_multiply(MTM,MT);
}

// Returns a matrix with is `M` with `V` appended as the last row (ie bottom)
function matrix_augment_row(M,V){

	// Make sure the lengths are the same (error)
	if( M[0].length !== V.length ){ return [[]] ;}

	// Copy the original matrix
	var R = [];
	var h = M.length;
	var i; var j; var w; var r;
	for(i=0; i<h; i+=1){
		R[i]=[];
		r = M[i]; //cache row
		w = r.length;
		for(j=0; j<w; j+=1){
			R[i][j] = r[j];
		}
	}
	
	// Append the vector as the final row
	//i += 1; //don't inc i
	R[i] = [];
	w = V.length;
	for(j=0; j<w; j+=1){
		R[i][j] = V[j];
	}
	
	return R;
}

// Returns a transposed copy of matrix `M`. (ie rotate M by 90* clockwise)
function matrix_transpose(M){
	
	// Create blank matrix of the appropraite size
	var R = [];
	var h = M.length;
	var w = M[0].length;
	var i=0, j=0;
	for(i=0; i<w; i+=1){ R[i] = []; }
	
	// Copy the elements into the transpose
	for(i=0; i<h; i+=1){
		for(j=0; j<w; j+=1){
			R[j][i] = M[i][j];
		}
	}
	
	// Return the transpose
	return R;
}

// Returns `ML + s*MR`. `ML` is the left matrix, `s` is a scalar that gets
// multiplied into `MR` the right matrix. You can use `s=1` for regular
// addition or `s=-1` for subtraction.
function matrix_add(ML,MR,s){
	
	// Make sure the heights agree (error)
	var h = ML.length;
	if(h!==MR.length){return [[]];}
	
	// Make sure the widths agree (error)
	var w = ML[0].length;
	if(w!==MR[0].length){return [[]];}
	
	// Perform the computation for each element
	var MO = [];
	var i=0, j=0;
	for(i=0; i<h; i+=1){
		MO[i] = [];
		for(j=0; j<w; j+=1){
			//console.info(ML[i][j] +"+"+ s +"*"+ MR[i][j] +"="+ (ML[i][j] + s*MR[i][j]) );
			MO[i][j] = ML[i][j] + s*MR[i][j];
		}
	}
	return MO;
}

// Returns `ML*MR` (matrix multiplication). `ML` is the matrix on the left,
// `MR` is tha matrix on the right.
function matrix_multiply(ML,MR){
	
	// Get the sizes (w=width, h=height) of the (l=left, r=right) matrix
	var lW = ML[0].length || 0;
	var lH = ML.length || 0;
	var rW = MR[0].length || 0;
	var rH = MR.length || 0;
	
	// If the dimensions don't agree: exit
	// The width of the left must equal the height of the right
	// because we multiply a row of the left with a col of the r
	if(lW !== rH){return [[]];}
	
	// The output matrix MO has lH rows and rW cols
	var MO = [];
	var i=0, j=0, ii=0, c=0;
	
	//iterate over each cell in MO
	for(i=0; i<lH; i++){
		//create the row
		MO[MO.length] = [];
		//compute the value of the cell
		for(j=0; j<rW; j++){
			c = 0;
			//multiply the row of ML with col of MR
			for(ii=0; ii<lW; ii++){
				//add values in c
				c += ML[i][ii]*MR[ii][j];
			}
			MO[i][j] = c;
		}
	}
	return MO;
}

// Returns the inverse of matrix `M`.
function matrix_inverse(M){
	// I use Guassian Elimination to calculate the inverse:
	// (1) 'augment' the matrix (left) by the identity (on the right)
	// (2) Turn the matrix on the left into the identity by elemetry row ops
	// (3) The matrix on the right is the inverse (was the identity matrix)
	// There are 3 elemtary row ops: (I combine b and c in my code)
	// (a) Swap 2 rows
	// (b) Multiply a row by a scalar
	// (c) Add 2 rows
	
	//if the matrix isn't square: exit (error)
	if(M.length !== M[0].length){return [[]];}
	
	//create the identity matrix (I), and a copy (C) of the original
	var i=0, ii=0, j=0, dim=M.length, e=0;
	var I = [], C = [];
	for(i=0; i<dim; i+=1){
		// Create the row
		I[I.length]=[];
		C[C.length]=[];
		for(j=0; j<dim; j+=1){
			
			//if we're on the diagonal, put a 1 (for identity)
			if(i===j){ I[i][j] = 1; }
			else{     I[i][j] = 0; }
			
			// Also, make the copy of the original
			C[i][j] = M[i][j];
		}
	}
	
	// Perform elementary row operations 
	for(i=0; i<dim; i+=1){
		// get the element e on the diagonal
		e = C[i][i];
		
		// if we have a 0 on the diagonal (we'll need to swap with a lower row)
		if(e===0){
			//look through every row below the i'th row
			for(ii=i+1; ii<dim; ii+=1){
				//if the ii'th row has a non-0 in the i'th col
				if(C[ii][i] !== 0){
					//it would make the diagonal have a non-0 so swap it
					for(j=0; j<dim; j++){
						e        = C[i][j]; //temp store i'th row
						C[i][j]  = C[ii][j];//replace i'th row by ii'th
						C[ii][j] = e;       //repace ii'th by temp
						e        = I[i][j]; //temp store i'th row
						I[i][j]  = I[ii][j];//replace i'th row by ii'th
						I[ii][j] = e;       //repace ii'th by temp
					}
					//don't bother checking other rows since we've swapped
					break;
				}
			}
			//get the new diagonal
			e = C[i][i];
			//if it's still 0, not invertable (error)
			if(e===0){return [[]];}
		}
		
		// Scale this row down by e (so we have a 1 on the diagonal)
		for(j=0; j<dim; j++){
			C[i][j] = C[i][j]/e; //apply to original matrix
			I[i][j] = I[i][j]/e; //apply to identity
		}
		
		// Subtract this row (scaled appropriately for each row) from ALL of
		// the other rows so that there will be 0's in this column in the
		// rows above and below this one
		for(ii=0; ii<dim; ii++){
			// Only apply to other rows (we want a 1 on the diagonal)
			if(ii===i){continue;}
			
			// We want to change this element to 0
			e = C[ii][i];
			
			// Subtract (the row above(or below) scaled by e) from (the
			// current row) but start at the i'th column and assume all the
			// stuff left of diagonal is 0 (which it should be if we made this
			// algorithm correctly)
			for(j=0; j<dim; j++){
				C[ii][j] -= e*C[i][j]; //apply to original matrix
				I[ii][j] -= e*I[i][j]; //apply to identity
			}
		}
	}
	
	//we've done all operations, C should be the identity
	//matrix I should be the inverse:
	return I;
}

// Returns the determinant (scalar) of square matrix `M`.
function matrix_determinant(M){
	// Use Guassian Elimination to turn the matrix into a triangular matrix
	// since the lower triangle is all 0's the determinant is just the
	// product of all of the diagonal elements times +1 or -1.
	// I forgot where I read about this but it seems to work.
	
	// If the matrix isn't square: exit (error)
	if(M.length !== M[0].length){ return; }
	// Start with a determinate of 1
	var det = 1;
	
	//create a copy of the original
	var i=0, ii=0, j=0, dim=M.length, e=0, C=[], f;
	for(i=0; i<dim; i++){
		C[C.length]=[];
		for(j=0; j<dim; j++){
			C[i][j] = M[i][j];
		}
	}
	
	// Perform elementary row operations 
	for(i=0; i<dim; i+=1){
		//get the element e on the diagonal
		e = C[i][i];
		
		// if we have a 0 on the diagonal (we'll need to swap with a lower row)
		if(e===0){
			//look through every row below the i'th row
			for(ii=i+1; ii<dim; ii+=1){
				//if the ii'th row has a non-0 in the i'th col
				if(C[ii][i] !== 0){
					//it would make the diagonal have a non-0 so swap it
					for(j=0; j<dim; j++){
						f        = C[i][j]; //temp store i'th row
						C[i][j]  = C[ii][j];//replace i'th row by ii'th
						C[ii][j] = f;       //repace ii'th by temp
					}
					//every swap requires us to the change the sign of det
					det *= -1;
					//don't bother checking other rows since we've swapped
					break;
				}
			}
			//get the new diagonal
			e = C[i][i];
			//if it's still 0, we have 0 det
			if(e===0){return 0;}
		}
		
		// We only care about the rows below this one 
		for(ii=i+1; ii<dim; ii+=1){
			// `e` is the value on the diagonal above this (ie same column)
			// Scale this value `C[ii][i]` by `e` to get the factor `f` that
			// can be used to make the element 0.
			f = C[ii][i] / e;
			
			// Subtract (the row above scaled by f) from (the current row)
			for(j=0; j<dim; j++){
				//chage the ii'th row by the i'th row times factor f
				C[ii][j] -= f*C[i][j];
			}
		}
		//we've got to multiply the diagonal elements
		//might as well do it as we calc them
		det *= e;
	}
	return det;
}

// Returns the `c`th column of matrix `M` as a column matrix
function matrix_get_col(M,c){
	var i=0, l=M.length;
	var a = [];
	for(i=0; i<l; i++){
		a[i] = [M[i][c]];
	}
	return a;
}

// Returns the `r`th row of matrix `M` as a row matrix
function matrix_get_row(M,r){
	var i=0, l=M[r].length;
	var a = [];
	a[0] = [];
	for(i=0; i<l; i++){
		a[0][i] = M[r][i];
	}
	return a;
}

// Returns a column matrix from an array of elements
function matrix_col_from_array(arr){
	var i=0, l=arr.length;
	var M = [];
	for(i=0; i<l; i++){
		M[i] = [arr[i]];
	}
	return M;
}

// Returns a column matrix from a vector
function matrix_col_from_vector(V){return matrix_col_from_array(V);}

// Returns an array created from the column matrix `M`
function matrix_col_to_array(M){
	var i=0, l=M.length;
	var a = [];
	for(i=0; i<l; i++){
		a[i] = M[i][0];
	}
	return a;
}
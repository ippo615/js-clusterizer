// Returns (`V` with scalar `s` inserted as the first element)
function vector_insert_first(V,s){
	var R = [s]; var l = V.length;
	for(var i=0; i<l; i+=1){
		R[R.length] = V[i];
	}
	return R;
}

// Returns (`V` with scalar `s` appended as the last element)
function vector_insert_last(V,s){
	var R = []; var l = V.length;
	for(var i=0; i<l; i+=1){
		R[R.length] = V[i];
	}
	R[R.length] = s;
	return R;
}

// Returns a copy of `V*s` (vector `V` scaled by scalar `s`)
function vector_scale(V,s){
	var l = V.length; var i;
	var R = [];
	for(i=0; i<l; i+=1){
		R[i] = V[i]*s;
	}
	return R;
}

// Returns the squared magnitude (scalar) of a vector. The squared magnitude
// is the sum of ( each element of the vector, squared ).
function vector_magnitude_squared(V){
	var d = V.length;
	var e; var sum=0; var i;
	
	// Compute the sum of (each element squared)
	for(i=0; i<d; i+=1){
		e = V[i];
		sum += e*e;
	}
	
	// return the squared magnitude
	return sum;
}

// Returns the magnitude (scalar) of a vector. This is also
// consided the length of the vector.
function vector_magnitude(V){
	return Math.sqrt(vector_magnitude_squared(V));
}
var vector_length = vector_magnitude;

// Returns the vector that is the result of applying unary operator (single
// argument function) `f` to each element of `V`. `f` is a function that should
// accept 1 argument (a number) and return a number.
function vector_unary(V,f){
	var R = [];
	var l = V.length;
	for(var i=0; i<l; i+=1){
		R[R.length] = f(V[i]);
	}
	return R;
}

// Returns VA+s*VB, (VA and VB are vectors, s is a scalar)
function vector_add(VA,VB,s){

	//if the vectors are different lengths...what do we do? (error)
	var al = VA.length; var bl = VB.length;
	if(al !== bl){ return; }
	
	//perform the computation
	var i; var r = [];
	for(i=0; i<al; i+=1){
		r[i] = VA[i] + s*VB[i];
	}
	
	// return the resultant vector
	return r;
}

// Returns a vector that represents the distance between vectors V1 and V2
// ie M1 = [4,5,6]; M2 = [1,2,3];
//    returns SQRT( (1-4)^2 + (5-2)^2 + (6-3)^2 );
function vector_distance(V1,V2){
	// Make sure the lengths are the same (error)
	if( V1.length !== V2.length ){ return ; }
	
	// Sum the squared differences between corresponding elements
	var i=0, l=V1.length, sum=0, d=0;
	for(i=0; i<l; i+=1){
		d = V2[i]-V1[i];
		sum += d*d;
	}

	// The distance is the square-root of the sum
	return Math.sqrt(sum);
}

// Return the product (a scalar) of A and B (ie `A dot B`)
function vector_dot(A,B){
	// make sure the dimensions are the same (error)
	var d = A.length;
	if( d !== B.length ){ return; }
	
	// iterate over each dimesion (element) in the vector add the product
	// of the corresponding elements in A and B to the sum
	var sum=0; var i;
	for(i=0; i<d; i+=1){
		sum += A[i]*B[i];
	}
	
	//return the sum
	return sum;
}

// Returns 1 if each element in VA equals the corresponding element in VB;
// otherwise, returns 0
function vector_is_equal(VA,VB){

	// if the vectors are different lengths then they're not equal
	var al = VA.length; var bl = VB.length;
	if(al !== bl){ return 0; }
	var i;
	
	// if any of the elements are different:
	// then the vectors can't be the same
	for(i=0; i<al; i+=1){
		if(VA[i] !== VB[i]){
			return 0;
		}
	}
	
	// we made it this far, the vectors are identical
	return 1;
}

// Returns the `r`th row of matrix `M` as a vector
function vector_from_matrix_row(M,r){
	var i=0; var l=M[r].length;
	var V = [];
	for(i=0; i<l; i++){
		V[i] = M[r][i];
	}
	return V;
	
}

// Returns the `c`th col of matrix `M` as a vector
function vector_from_matrix_col(M,c){
	var V = [];	
	var l = M.length;
	for(var i=0; i<l; i+=1){
		V[V.length] = M[i][c];
	}
	return V;
}

// Returns a row matrix of the elements of vector `V`
function vector_to_row_matrix(V){
	var R = []; 
	R[0]=[];
	var i; var l = V.length;
	for(i=0; i<l; i+=1){
		R[0][i] = V[i];
	}
	return R;
}
/************************************************************************[LDF]*
 *  ldf_calculate_perceptron(ClassA,ClassB,MaxIter)
 *  ClassA is an array of cluster objects
 *  ClassB is an array of cluster objects
 *  ClassB will be negated
 *  MaxIter is the maximum number of iteratrions to allow
 */
function ldf_calculate_perceptron(ClassA,ClassB,MaxIter){

	// Determine sized and number of points
	var d = cluster_get_dimension(ClassA);
	var j, i, la, lb;
	la = cluster_get_num_points(ClassA);
	lb = cluster_get_num_points(ClassB);

	// Set the perceptron P and previous perceptron PP values to 0;
	var P = []; var PP = [];
	for(i=0; i<d; i+=1){
		P[i]  = 0; 
		PP[i] = 0;
	}

	//only loop for the maximum number of iterations
	for(i=0; i<MaxIter; i++){
		
		// Store the previous perceptron values
		PP = P;

		// Find the elements in A that were classified incorrectly
		for(j=0; j<la; j++){

			// Adjust the perceptron if a point is wrongly classified
			// I include the 0 case too; otherwise when you start
			// with 0,0 nothing would happen
			if(vector_dot(PP,cluster_get_point(ClassA,j)) <= 0){

				// You could change 1 to some 'weight' which
				// controls how fast/slow the perceptron learns
				P = vector_add(P,cluster_get_point(ClassA,j),1);

			}
		}

		// Find the elements in B that were classified incorrectly
		for(j=0; j<lb; j++){

			//instead of negating each element of ClassB and
			//checking if the dot product is negative:
			//just check if this dot product is positive
			if(vector_dot(PP,cluster_get_point(ClassB,j)) >= 0){

				//we need to subtract these values from P
				P = vector_add(P,cluster_get_point(ClassB,j),-1);
			}
		}

		// If the perceptron is the same, it works!
		if(vector_is_equal(P,PP)){break;}
	}

	// The actual 'dividing line' is perpendicular to the perceptron;
	// however, determining the line in n-dimensions is hard.
	return P;
}

/************************************************************************[LDF]*
 *  ldf_calculate_perceptronAug(ClassA,ClassB,Opt)
 *  ClassA is an array of cluster objects
 *  ClassB is an array of cluster objects
 *  ClassB will be negated
 *  Opt is 'options' is an object: has maxIter propterty
 */
function ldf_calculate_perceptronAug(ClassA,ClassB,Opt){

	// Determine sized and number of points
	// Note this is 'agumented' so add an extra dimension
	var d = cluster_get_dimension(ClassA)+1;
	var j, i, la, lb,e ;
	la = cluster_get_num_points(ClassA);
	lb = cluster_get_num_points(ClassB);

	// Set the perceptron P and previous perceptron PP values to 0;
	var P = []; var PP = [];
	for(i=0; i<d; i+=1){
		P[i]  = 0; 
		PP[i] = 0;
	}

	//only loop for the maximum number of iterations
	var MaxIter = Opt.maxIter || 100;
	for(i=0; i<MaxIter; i++){

		//store the previous perceptron values
		PP = P;

		// Adjust the perceptron based on the misclassified elements in A
		for(j=0; j<la; j++){

			// Remeber to augment the data point and correct the perceptron
			e = vector_insert_first(cluster_get_point(ClassA,j),1);
			if(vector_dot(PP,e) <= 0){
				P = vector_add(P,e,1);
			}
		}

		// Adjust the perceptron based on the misclassified elements in A
		for(j=0; j<lb; j++){

			// Remeber to augment the data point and correct the perceptron
			e = vector_insert_first(cluster_get_point(ClassB,j),1);
			if(vector_dot(PP,e) >= 0){
				P = vector_add(P,e,-1);
			}
		}

		// If the perceptron is the same, it works!
		if(vector_is_equal(P,PP)){break;}
	}

	// The actual 'dividing line' is perpendicular to the perceptron;
	// however, determining the line in n-dimensions is hard.
	return P;

}

/************************************************************************[LDF]*
 *  MSE_LDF(ClassA,ClassB,opt)
 *  ClassA is an array of cluster objects
 *  ClassB is an array of cluster objects
 *  Opt is a options object:
 *    Opt.isFisher == 1 -> fisher ldf
 *    Opt.isFisher == 0 -> normal ldf
 */
function ldf_calculate_fisher(ClassA,ClassB,Opt){Opt.isFisher=1; return _ldf_calculate_ldf_mse(ClassA,ClassB,Opt);}
function ldf_calculate_mse(ClassA,ClassB,Opt){Opt.isFisher=0; return _ldf_calculate_ldf_mse(ClassA,ClassB,Opt);}
function _ldf_calculate_ldf_mse(ClassA,ClassB,Opt){

	// We're solving: Y*A=B -> A=B/Y (note: A,B are vectors, Y is matrix)

	// determine the number of elements in each class
	var la = cluster_get_num_points(ClassA);
	var lb = cluster_get_num_points(ClassB);

	// Create the Y matrix, dim wide by la+lb tall
	// Augment all of the data points by 1; then negate class B
	// Then insert each element as a row in Y (A on top, B on bottom)
	var Y = [];
	var i, j, e;
	for(i=0; i<la; i+=1){
		e = cluster_get_point(ClassA,i);
		Y[i] = vector_insert_first(e,1);
	}
	for(j=0; j<lb; j+=1, i+=1){
		e = cluster_get_point(ClassB,j);
		Y[i] = vector_scale(vector_insert_first(e,1),-1);
	}

	// Create the b vector, Note the b vector is slighly different
	// if this is the Fisher variant of the MSE
	var isFisher = Opt.isFisher || 0;
	var tmp = []; var n = la+lb;
	if( isFisher ){

		// For the first (# of elements in A)
		// the element of this vector is the ratio of:
		// (# of elements in A):(total # elements)
		for(i=0; i<la; i+=1){
			tmp = vector_insert_last(tmp,la/n);
		}

		// For the remaining (# of elements in B)
		// the element of this vector is the ratio of:
		// (# of elements in B):(total # elements)
		for(j=0; j<lb; j+=1){
			tmp = vector_insert_last(tmp,lb/n);
		}
	}else{

		// For the regular MSE, b is all 1's
		for(i=0; i<la+lb; i+=1){
			tmp = vector_insert_last(tmp,1);
		}
	}

	// Since Y is probably not be square, we need to compute the
	// psuedo inverse.
	var b = matrix_col_from_vector(tmp);
	var YPS = matrix_psuedo_inverse(Y);

	// Finally compute and return A
	var A = matrix_multiply(YPS,b);
	return vector_from_matrix_col(A,0);
}

/************************************************************************[LDF]*
 *  ldf_calculate_ho_kashyap(ClassA,ClassB,Opt)
 *  ClassA is an array of cluster objects
 *  ClassB is an array of cluster objects
 *  Opt
 *    maxIter -> maximum number of iterations to allow
 *    maxErr  -> maximum error to allow
 */
function ldf_calculate_ho_kashyap(ClassA,ClassB,Opt){
	
	var maxIter = Opt.maxIter || 999;
	var maxErr  = Opt.maxErr || 1.0;
	
	//create augmented versions of the classes
	var CAA = cluster_copy_augment(ClassA,1);
	var CBA = cluster_copy_augment(ClassB,1);
	
	//begin creating the Y matrix
	var Y = vector_to_row_matrix(cluster_get_point(CAA,0));
	
	// Add the remaining elements from class A
	// i=1 because we already added the 0 element above
	var i; var la = cluster_get_num_points(CAA);
	for(i=1; i<la; i+=1){
		Y = matrix_augment_row(Y,cluster_get_point(CAA,i));
	}
	
	// Add the negated elements from class B to the matrix
	var j; var lb = cluster_get_num_points(CBA); var pt;
	for(j=0; j<lb; j+=1){
		pt = cluster_get_point(CBA,j);
		Y = matrix_augment_row(Y,vector_scale(pt,-1));
	}

	//compute Ydagger (the psuedo inverse)
	var YD = matrix_psuedo_inverse(Y);


	//create the b vector with all 1's
	var tmp = [];
	var n = la+lb;
	for(i=0; i<n; i+=1){ tmp[i] = 1; }
	var b = matrix_col_from_array(tmp);

	// We only care about positive errors (this will filter those)
	var getPositive = function(a){
		if(a>0){return a;}else{return 0;}
	};
	
	var a; var e; var ep; var YA; var bV; var t;
	while(maxIter--){

		// Compute a (the answer)
		a = matrix_multiply(YD,b);

		// Compute the error, e
		YA = vector_from_matrix_col(matrix_multiply(Y,a),0);
		bV = vector_from_matrix_col(b,0);
		e  = vector_add(YA,bV,-1);

		// Get the positive part of the error
		ep = vector_unary(e,getPositive);

		// We should be minimizing:
		// f(t) = ||Y*a-(b-t*(Y*a-b))||^2
		// f(t) = ||e+te||^2 for t
		// then b = b-t*b -> where did i get this? class notes
		// or b = a+2t*ep -> www.cedar.buffalo.edu/~srihari/CSE555/Chap5.Part3.pdf
		// but I use b = b+2t*ep -> um... cuase it makes sense and can be computed
		// t appears to be a number
		// asssume t=1
		t = 0.3;
		bV = vector_add(bV,ep,2*t);
		b = matrix_col_from_array(bV);

		// if error is small enough, return a solution
		if( vector_magnitude_squared(e) < maxErr){
			return a;
		}
	}

	// even though the error is too large return our answer
	return a;
}

/**********************************************************************[SCORE]*
 * Calucates the sum of the scatter within matrices of a set of clusters
 */
function scatter_within_sum_matrix(clusters){

	// Create a blank scatter total matrix
	var dim = cluster_get_dimension(clusters[0]);
	var ST = matrix_make_zero(dim,dim);

	// Add each cluster's scatter within matrix to the total matrix
	var l = clusters.length; var i; var SW;
	for(i=0; i<l; i++){
		SW = scatter_within_make_matrix(clusters[i]);
		ST = matrix_add(ST,SW,1);
	}

	return ST;
}

/**********************************************************************[SCORE]*
 * Calculates the scatter between matrix for a set of clusters
 */
function scatter_between_make_matrix(clusters){

	// Create a blank scatter total matrix
	var dim = cluster_get_dimension(clusters[0]);
	var ST = matrix_make_zero(dim,dim); var SW;

	// Get a count of the total # of points
	// vecSum will serve as the center of the entire data set
	var totalPts=0, vecSum=[], i;
	for(i=0; i<dim; i++){vecSum[i]=0;}


	// Find the center of our entire data set and the total number of
	// points in our data set (ie pretend ALL clusters are 1 cluster)
	var l = clusters.length;
	var c; var j; var pts;
	var dV; var dVc; var dVr; var dVm;
	for(i=0; i<l; i++){
		c = clusters[i];
		pts = cluster_get_num_points(c);
		for(j=0; j<pts; j+=1){
			totalPts += 1;
			vecSum = vector_add(vecSum,cluster_get_point(c,j),1);
		}
	}

	// Divide the sum by the total # pts to get the center of all data
	vecSum = vector_scale(vecSum,1/totalPts);

	// Compute the scatter between matrix (like the scatter within)
	// Iterate over all clusters
	for(i=0; i<l; i++){
		c = clusters[i];
		SW = matrix_make_zero(dim,dim);

		// Iterate over every point in the cluster
		pts = cluster_get_num_points(c);
		for(j=0; j<pts; j+=1){

			// Compute the difference from the center
			dV = vector_add(vecSum,cluster_get_point(c,j),-1);

			// Convert to a col matrix and row matrix for mult
			dVc = matrix_col_from_array(dV);
			dVr = matrix_transpose(dVc);

			// Multiply column times row to get a matrix
			dVm = matrix_multiply(dVc,dVr);

			// Add it to the scatter matrix
			SW = matrix_add(SW,dVm,1);
		}

		// Add the scatter matrix for the cluster to the total 
		ST = matrix_add(ST,SW,1);
	}

	return ST;
}

/********************************************************************[SCATTER]*
 * Calucates the scatter within matrix for a cluster
 */
function scatter_within_make_matrix(cluster){

	// Start with a blank scatter matrix
	var dim = cluster_get_dimension(cluster);
	var SW = matrix_make_zero(dim,dim);

	// Get the mean of the cluster
	var m = cluster_calc_center(cluster);

	// The scatter within matrix is the sum of the difference between
	// each point and the mean but in matrix form. 
	var pts = cluster_get_num_points(cluster);
	var j; var dV; var dVc; var dVr; var dVm;
	for(j=0; j<pts; j+=1){

		// Compute the difference from the center
		dV = vector_add(m,cluster_get_point(cluster,j),-1);

		// Convert to a col matrix and row matrix for mult
		dVc = matrix_col_from_array(dV);
		dVr = matrix_transpose(dVc);

		// Multiply column times row to get a matrix
		dVm = matrix_multiply(dVc,dVr);

		// Add it to the scatter matrix
		SW = matrix_add(SW,dVm,1);
	}

	return SW;
}

/* XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX *
 * X X X X X X X X X X X X X X X X X X X X         Experimental Stuff Below  *
 * XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX */

/************************************************************************[LDF]*
 *  ldf_calculate_scatter(ClassA,ClassB,Opt) -- unfinished
 *  ClassA is an array of cluster objects
 *  ClassB is an array of cluster objects
 *  Opt does nothing currently
 */
function ldf_calculate_scatter(ClassA,ClassB,Opt){
	// Augment the classes (this creates a row of 0's in the scatter
	// matrix, making inversion impossible)
	var CAA = cluster_copy_augment(ClassA,1);
	var CBA = cluster_copy_augment(ClassB,1);

	// Create the scatter within matrices
	var SWA = scatter_within_make_matrix(CAA);
	var SWB = scatter_within_make_matrix(CBA);

	// The overall scatter within matrix
	var SW = matrix_add(SWA,SWB,1);

	// Compute the scatter between the augmented classes
	// var SB = scatter_between_make_matrix([CAA,CBA]);

	// Confusion!
	//console.info("SB"); //console.info(SB);
	//add the scatter between to the scatter within?????
	//according to my notes, I shouldn't add SB
	//SW = matrix_add(SW,SB,1);
	//console.info("SW"); //console.info(SW);

	// Invert the scatter matrix (impossible!)
	var SWI = matrix_inverse(SW);
	//console.info("SWI"); //console.info(SWI);

	// Compute the vector from B to A
	var mA = cluster_calc_center(CAA);
	var mB = cluster_calc_center(CBA);
	var wo = vector_add(mA,mB,-1);

	// Compute the w
	var w = matrix_multiply(SWI,matrix_col_from_array(wo));

	//convert to vector
	//proves the 'direction' of the line is correct
	//but the offset caused by augmenting is wrong
	//
	var r = vector_from_matrix_col(w,0);
	//r.splice(0,1,0);
	//remove the 'offset' from the line
	//ie make it pass through the origin
	r.splice(0,1,0);

	//compute the middle of the data set
	var m = vector_add(mA,mB,1);
	m = vector_scale(m,0.5);

	//the line is the dot prod of:
	//ie [1, xmid, ymid] [OFFSET, X, Y] = 0
	//therefore OFFSET = -[1,X,Y]*[0,xmid,ymid]
	//so put that computation in the offset
	r.splice(0,1,-vector_dot(r,m));

	return r;
}

/************************************************************************[LDF]*
 *  Support_Vector_Machine(ClassA,ClassB,MaxIter) -- unfinished
 */
function ldf_calculate_svm(ClassA,ClassB){

	// Augment the classes
	var CAA = cluster_copy_augment(ClassA,1);
	var CBA = cluster_copy_augment(ClassB,1);

	//polulate the Y matrix with the augment classes
	var la = cluster_get_num_points(CAA);
	var lb = cluster_get_num_points(CBA);
	var i=0; var Y = [];
	for(i=0; i<la; i+=1){Y = matrix_augment_row(Y,cluster_get_point(CAA,i));}
	for(i=0; i<lb; i+=1){Y = matrix_augment_row(Y,cluster_get_point(CBA,i));}

	// Populate Z with +1 for class A elements, -1 for class B
	var Z = [];
	i=0;
	for(i=0; i<la; i+=1){ Z[Z.length] =  1; }
	for(i=0; i<lb; i+=1){ Z[Z.length] = -1; }

	var n = la+lb-1;	

	//console.info(Y);
	//compute L:
	/*
	var j=0;
	var L="";
	for(i=0; i<n; i+=1){
		for(j=0; j<n; j+=1){
			L += "-0.5*l"+i+"*l"+j;
			L += "*"+Z[i]*Z[j];
			//console.info(matrix_get_row(Y,i));
			//console.info(matrix_get_row(Y,j));
			//console.info(matrix_multiply(matrix_get_row(Y,i),matrix_transpose(matrix_get_row(Y,j)),1));
			L += "*"+matrix_multiply(matrix_get_row(Y,i),matrix_transpose(matrix_get_row(Y,j)),1)[0][0];
			//(-1/2)*l[i]*l[j]*z[i]*z[j]*matrix_multiply(matrix_transpose(Y[i]),Y[j],1);
		}
	}
	*/
	
	// Instead of computing L we're going to compute an array of dL/di
	// L is "some number" times li PLUS "some (other?) number" times li*li PLUS some constants
	// dL/di is "sume number" PLUS "double some (other?) number" times li
	var j=0;
	var L=[];
	for(i=0; i<n; i+=1){
		L[i] = "";
		for(j=0; j<n; j+=1){
			//if we've got the li*li case:
			if(j===i){
				//the 2 cancels with the 0.5 yielding:
				L[i] += "-1*l"+j;
			}else{
				//some number (remember lj is assumed constant for given li)
				L[i] += "-0.5*l"+j;
			}
			//this "number" is always the same
			L[i] += "*"+Z[i]*Z[j];
			L[i] += "*"+matrix_multiply(matrix_get_row(Y,i),matrix_transpose(matrix_get_row(Y,j)),1)[0][0];
			//
			//(-1/2)*l[i]*l[j]*z[i]*z[j]*matrix_multiply(matrix_transpose(Y[i]),Y[j],1);
		}
	}

	// Create the "names" of the symbolic variables used in the equation
	var varNames = [];
	for(i=0; i<n; i+=1){ varNames[i] = "l"+i; }

	// Create a system of equations with the unknowns (varNames)
	// and the derivatives (L)
	var F = CreateSystemOfEquations(L,varNames);

	// The initial guess is [0]
	var guess = [];
	for(i=0; i<n; i+=1){ guess[i] = 10000; }

	// Solve the problem via newton's method
	var S = NewtonSolve(F,guess,0.1,0.001,250);

	//maximize L (or minimize -L) to get all of the l or lambdas

	return S;
}

/**************************************************************************************************
 * CreateSystemOfEquations(funcStrArr,variablesStrArr) -> Creates "row" matrix
 * Creates an array of functions that can be evaluated for any value of their arguments
 * funcStrArr is an array of strings that represents the functions
 * variablesStrArr is an array of strings that represents the variables
 * ie: F = CreateSystemOfEquations(["X+Y","X*Y"],["X","Y"]);
 *     F[0]([1,2]) == 1+2 == 3
 *     F[1]([1,2]) == 1*2 == 2
 */
function CreateSystemOfEquations(funcStrArr,variablesStrArr){
	var i=0, l=funcStrArr.length;
	var funcs = [];
	for(i=0; i<l; i++){
		funcs[funcs.length] = CreateNumericFunction(funcStrArr[i],variablesStrArr);
	}
	return funcs;
}
/****************************************************************************************[NUM LIB]*
 * The function creates a function that can be evaluated for any value of it's arguments
 * funcString is a string that represents the functions
 * varStrArr is an array of strings that represent the variables
 * ex: F1 = CreateNumericFunction("X+Y",["X","Y"]);
 *     F1([2,3]) == 2+3 == 5;
 * ex: G1 = CreateNumericFunction("DOG+CAT*DOG",["DOG","CAT"]);
 *     G1([2,5]) == 2+5*2 == 12;
 */
function CreateNumericFunction(funcString,varStrArr){

	// Create a string to be "eval"ed that will create a function called 'x'
	// ie "var x = function($){return $[0]+$[1]*$[2];};"
	var evalFunc = "x = function($){";

	// We want to return the value of the function expression
	evalFunc += " return ";

	// Replace all of the variable names with array indexes
	// ie X->$[0], Y->$[1]
	var i=0, l=varStrArr.length;
	var tmp = funcString;

	// TEMP HACK:
	// iterate backwards because when we generate the l's
	// in numberic order ie l1, l2, l3, so if we replace in that order
	// l1->$[1] will also replace l12->$[1]2 which is a syntax error
	for(i=l; i>=0; i--){
		while(tmp.indexOf(varStrArr[i]) > -1){
			tmp = tmp.replace(varStrArr[i],"$["+i+"]");
		}
	}
	evalFunc += tmp;
	evalFunc += ";};";

	// Create the function by evaling, then return it:
	var x;
	eval(evalFunc);
	return x;
}

function CalcJacobian(funcArr,guessVec,delta){

	var J=[];
	var F=funcArr;
	var fl=F.length;
	var vl=guessVec.length;
	var i=0, j=0, k=0;
	var Xr, Xl;

	// For each function
	for(i=0; i<fl; i++){
		J[i] = [];

		// For each dimension in the guess vector
		for(j=0; j<vl; j++){
			//J[i,j] = dF[i]/dx[j];

			// Create a vector of the point we are testing
			Xr = matrix_col_to_array(guessVec);
			Xl = matrix_col_to_array(guessVec);

			// Move a component just to the right and left
			Xr[j] += delta;
			Xl[j] -= delta;

			// Numerically approximate the dervitive
			J[i][j] = ( F[i](Xr)-F[i](Xl) )/(2*delta);
		}
	}

	return J;
}
/**************************************************************************************************
 *
 */
function ui_decimal_to_fraction(decimal,maxDenom){
	var numerator   = Math.floor(decimal*maxDenom);
	var denominator = Math.floor(maxDenom);
	
	// Compute the greatest common denominator (GCD)
	var gcd = denominator;
	var num = numerator;
	var temp = gcd;
	while(num !== 0){
		temp = num;
		num = gcd % num;
		gcd = temp;
	}

	// Return the numerator & denominaotr scaled by GCD
	return {N: numerator/gcd, D: denominator/gcd};
}

/**************************************************************************************************
 * xp = x - (J(x)^-1)*f(x)
 */
function NewtonSolve(SYS,Guess,Delta,MaxError,MaxIter){

	// Calcute the Jacobian for the 1st guess
	var Jo = CalcJacobian(SYS,matrix_col_from_array(Guess),Delta);

	//console.info("Initial Jacobian:");
	//console.info(Jo);
	//alert(Jo);

	// Calcute the inverse Jacobian
	var Ji = matrix_inverse(Jo);

	//console.info("Initial Inverse Jacobian:");
	//console.info(Ji);
	//alert(Ji);

	// Compute the F(x) vector
	var F = []; var i=0; var l=SYS.length;
	for(i=0; i<l; i++){
		F[i] = SYS[i](Guess);
	}
	//console.info("Initial F array:");
	//console.info(F);

	// Compute: (J(x)^-1)*f(x)
	var JiF = matrix_multiply(Ji,matrix_col_from_array(F));
	//console.info("Inverse jacobian times F:");
	//console.info(JiF);
	//alert(JiF);

	// Compute: x - (J(x)^-1)*f(x)
	var XP = matrix_add(matrix_col_from_array(Guess),JiF,-1);
	//console.info("x - (J(x)^-1)*f(x)");
	//console.info(XP);


	Jo = CalcJacobian(SYS,XP,Delta);
	Ji = matrix_inverse(Jo);
	F = [];
	for(i=0; i<l; i++){F[i] = SYS[i](matrix_col_to_array(XP));}
	JiF = matrix_multiply(Ji,matrix_col_from_array(F));
	var X = matrix_add(XP,JiF,-1);

	// Set a maximum number of iterations
	var iter = 0;
	var tE = MaxError;
	var error = MaxError;

	// if the difference between X_n and X_n-1 is > than the precision
	// and we haven't exceded the number of iterations
	while(error >= MaxError && iter++ < MaxIter){
		// Store the old X
		XP = X;

		// Compute the new X
		Jo = CalcJacobian(SYS,XP,Delta);
		Ji = matrix_inverse(Jo);
		for(i=0; i<l; i++){F[i] = SYS[i](matrix_col_to_array(XP));}
		JiF = matrix_multiply(Ji,matrix_col_from_array(F));
		X = matrix_add(XP,JiF,-1);

		// Sum the absolute value of each equation's error
		error = 0;
		for(i=0; i<l; i++){
			error += Math.abs(SYS[i](matrix_col_to_array(X)));
		}
	}

	// Convert to fractions (to double check)
	var round = [];
	var md = 50;
	for(i=0; i<l; i++){
		round[i] = ui_decimal_to_fraction(X[i][0],md);
	}

	return {
		SOLUTION: X,
		SOL_FRAC: round,
		ITERATIONS: iter,
		ERROR: error
	};
}
/**********************************************************************[SCORE]*
 * Calucates the SAVGE (Sum of Average Squared Errors) of an array of clusters
 */
function score_avg_error(clusters){

	// Start with 0 sum of squared errors
	var SSE = 0;

	// c will be a cluster, m is the middle of cluster
	var c; var m; var d;
	var i; var j; var pts; var l;
	l = clusters.length;

	//iterate over every cluster
	for(i=0; i<l; i+=1){

		c = clusters[i];
		m = cluster_calc_center(c);

		//iterate over every point in the cluster
		pts = cluster_get_num_points(c);
		for(j=0; j<pts; j+=1){

			// Difference between the middle and this point
			d = vector_add(m,cluster_get_point(c,j),-1);

			// Should I divide by the number of points in the
			// cluster? This would allow clusters with varying
			// number of points to have the same influence
			// SSE += vector_magnitude_squared(d)/pts;
			SSE += vector_magnitude_squared(d);

		}
	}

	// Should I divide by the total number of clusters?
	// return SSE/l;
	return SSE;
}
/**********************************************************************[SCORE]*
 * Calucates the SMAXE (Sum of MAXimum squared Errors) of an array of clusters
 */
function score_max_error(clusters){

	var SSE = 0;
	var c; var m; var dV; var dM; var max;
	var i; var j; var pts; var l;
	l = clusters.length;

	// Add the maximum error from each cluster
	for(i=0; i<l; i+=1){

		c = clusters[i];
		m = cluster_calc_center(c);

		// Start the max with a valid point
		dV = vector_add(m,cluster_get_point(c,0),-1);
		max = vector_magnitude_squared(dV);

		// Find the maximum error for this cluster
		pts = cluster_get_num_points(c);
		for(j=0; j<pts; j+=1){

			// If we have a new max, remember it
			dV = vector_add(m,cluster_get_point(c,j),-1);
			dM = vector_magnitude_squared(dV);
			if(dM > max){
				max = dM;
			}
		}

		// Add the max error
		SSE += max;
	}

	// Divide the sum by the total number of clusters?
	// return SSE/l;
	return SSE;
}
/**********************************************************************[SCORE]*
 * Calucates the SMINE (Sum of MINimum squared Errors) of an array of clusters
 */
function score_min_error(clusters){

	var SSE = 0;
	var c; var m; var dV; var dM; var min;
	var i; var j; var pts; var l;
	l = clusters.length;

	// Add the minimum error from each cluster
	for(i=0; i<l; i+=1){

		c = clusters[i];
		m = cluster_calc_center(c);

		dV = vector_add(m,cluster_get_point(c,0),-1);
		min = vector_magnitude_squared(dV);

		// Find the minimum error in the cluster
		pts = cluster_get_num_points(c);
		for(j=0; j<pts; j+=1){

			// If we have a new min, remember it
			dV = vector_add(m,cluster_get_point(c,j),-1);
			dM = vector_magnitude_squared(dV);
			if(dM < min){
				min = dM;
			}

		}

		// Update the sum of the errors
		SSE += min;
	}

	// Divide the sum by the total number of clusters?
	// return SSE/l;
	return SSE;
}
/**********************************************************************[SCORE]*
 * Calucates the SMEDE (Sum of MEDian Squared Errors) of an array of clusters
 */
function score_med_error(clusters){

	var SSE = 0;
	var c; var m; var d; var errors = [];
	var i; var j; var pts; var l;
	l = clusters.length;

	// Add the median error for each cluster
	for(i=0; i<l; i+=1){

		// start with 0 errors for this cluster
		errors = [];
		c = clusters[i];
		m = cluster_calc_center(c);

		// Create a list of all the errors in this cluster
		pts = cluster_get_num_points(c);
		for(j=0; j<pts; j+=1){
			d = vector_add(m,cluster_get_point(c,j),-1);
			errors[errors.length] = vector_magnitude_squared(d);
		}

		// Add the median error to the SSE
		SSE += stat_median(errors);
	}

	// Divide by the total number of clusters?
	//return SSE/l;
	return SSE;
}
/**********************************************************************[SCORE]*
 * Calucates the trace of the scatter matrix which is the sum
 * of the scatter within matrices and scatter between matrix
 */
function score_trace_scatter_matrix(clusters){
	var SW = scatter_within_sum_matrix(clusters);
	var SB = scatter_between_make_matrix(clusters);
	var ST = matrix_add(SW,SB,1);
	//return matrix_trace(ST)/(ST.length*ST.length);
	return matrix_trace(ST);
}
/**********************************************************************[SCORE]*
 * Calucates the trace of the scatter matrix which is the sum
 * of the scatter within matrices and scatter between matrix
 */
function score_trace_sw_sb(clusters){
	var SW = scatter_within_sum_matrix(clusters);
	var SB = scatter_between_make_matrix(clusters);
	var ST = matrix_multiply(matrix_transpose(SW),SB,1);
	//return matrix_trace(ST)/(ST.length*ST.length);
	return matrix_trace(ST);
}
function score_trace_sw_st(clusters){
	var SW = scatter_within_sum_matrix(clusters);
	var SB = scatter_between_make_matrix(clusters);
	var ST = matrix_add(SW,SB,1);
	var R  = matrix_multiply(SW,matrix_transpose(ST),1);
	//return matrix_trace(R)/(R.length*R.length);
	return matrix_trace(R);
}
function score_det_sw(clusters){
	var SW = scatter_within_sum_matrix(clusters);
	return matrix_determinant(SW);
}
function score_det_swSt(clusters){
	var SW = scatter_within_sum_matrix(clusters);
	var SB = scatter_between_make_matrix(clusters);
	var ST = matrix_add(SW,SB,1);
	var R = matrix_multiply(matrix_transpose(ST),SW,1);
	return matrix_determinant(R);
}

/**********************************************************************[STATS]*
 * Returns the median of an array
 * data is an array of numbers
 */
function stat_median(data){
	//determine the index of the middle element
	var half = data.length/2;
	//if we've only got 1 point, return that point;
	if(half === 0.5){return data[0];}
	//round the halfway point down
	var mid = Math.floor(half);
	//sort the array is ascending order (WARNING: ANNIHILATES DATA)
	//so I'm being lazy, TODO: implement quicksort or shell or insertion
	data.sort(function(x,y){return x-y;});
	//if the rounded middle is the same as the actual middle
	if (mid === half){
		//we've got a middle data point
		return 0+data[half];
	}else{
		//return the average of the middles
		return 0+(data[mid] + data[mid+1])/2;
	}
}
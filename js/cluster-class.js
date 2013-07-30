/********************************************************************[CLUSTER]*
 * Returns a new cluster object with the following properties
 * CENTROID: array of numbers (ie row vector) describing the center of mass
 * SD: standard deviation vector (ie [x_sd, y_sd, ...])
 * VAR: variance vector (ie standard deviation vector squared)
 * POINTS: [[x1,y1,...],[x2,y2,...]] - collection of data points in cluster
 * DIM: integer describing the dimensionality of the data
 * DIRTY: has the data changed and the stats not been updated? yes=1
 */
function cluster_create(DIM){
	return {
		CENTROID:[],
		SD:[],
		VAR:[],
		POINTS:[],
		DIM:DIM,
		DIRTY:1
	};
}

/********************************************************************[CLUSTER]*
 * Returns a cluster object that is a copy of `C`
 * It increases the dimension by 1 and sets the new fist element
 * of each point to `s`
 */
function cluster_copy_augment(C,s){
	var NC = cluster_create(cluster_get_dimension(C)+1);
	var l = cluster_get_num_points(C);
	var pt; var i;
	for(i=0; i<l; i+=1){
		pt = cluster_get_point(C,i);
		cluster_add_point(NC,vector_insert_first(pt,s));
	}
	return NC;
}

/********************************************************************[CLUSTER]*
 * Returns a copy of cluster object `C`
 * It also scales all of the points in the copy by `s`
 */
function cluster_copy_scale(C,s){
	var NC = cluster_create(cluster_get_dimension(C));
	var l = cluster_get_num_points(C);
	var pt; var i;
	for(i=0; i<l; i+=1){
		pt = cluster_get_point(C,i);
		cluster_add_point(NC,vector_scale(pt,s));
	}
	return NC;
}

/********************************************************************[CLUSTER]*
 * Adds a point to the clusters array (no return value)
 * CLUSTER is the cluster that the point should be added to
 * POINT is an array of numbers representing the data point to be added
 */
function cluster_add_point(CLUSTER,POINT){
	// Make sure the POINT has the same dimensionality as the CLUSTER
	if(POINT.length === CLUSTER.DIM){
		CLUSTER.POINTS[CLUSTER.POINTS.length] = POINT;
		CLUSTER.DIRTY = 1;
	}
}

/********************************************************************[CLUSTER]*
 * Returns the INDEX'th (0-based) data point in CLUSTER
 */
function cluster_get_point(CLUSTER,INDEX){
	if(CLUSTER.POINTS.length < INDEX){return [];}
	return CLUSTER.POINTS[INDEX];
}

/********************************************************************[CLUSTER]*
 * Returns the number of data points in this cluster
 */
function cluster_get_num_points(CLUSTER){
	return CLUSTER.POINTS.length;
}

/********************************************************************[CLUSTER]*
 * Returns the dimension of the data points in this cluster
 */
function cluster_get_dimension(CLUSTER){
	return CLUSTER.DIM;
}

/********************************************************************[CLUSTER]*
 * Removes all the points from the cluster
 */
function cluster_clear_points(CLUSTER){
	CLUSTER.POINTS = [];
	CLUSTER.DIRTY  = 1;
}

/********************************************************************[CLUSTER]*
 * Computes the CLUSTER's center and returns the center
 * CLUSTER is the cluster that you want to know the center of
 */
function cluster_calc_center(CLUSTER){
	
	// The centroid will start as 0
	var mid = [];
	var i; var l = CLUSTER.DIM;
	for(i=0; i<l; i++){mid[i] = 0;}

	// Compute the average of all of dimensions of all the points	
	var PTS = CLUSTER.POINTS;
	l = PTS.length;
	var scale = 1/l;
	for(i=0; i<l; i+=1){
		mid = vector_add(mid,PTS[i],scale);
	}
	
	// Update and return the data
	CLUSTER.CENTROID = mid;
	CLUSTER.DIRTY = 0;
	return CLUSTER.CENTROID;
}

/********************************************************************[CLUSTER]*
 * Returns CLUSTER's variance (vector)
 * CLUSTER is the cluster that you want to know variance
 */
function cluster_calc_variance(CLUSTER){

	// Get the mean of each dimension
	var m = cluster_calc_center(CLUSTER);

	// Calculate the variance for each dimension
	var v = []; var i=0; var j; var diff;
	var d   = cluster_get_dimension(CLUSTER);
	var pts = cluster_get_num_points(CLUSTER);
	for(i=0; i<d; i+=1){

		// The variance is the average( ( difference from the mean) squared )
		// ie: sum( (x[i]-mean)^2 ) / n
		v[i] = 0;
		for(j=0; j<pts; j+=1){
			diff = cluster_get_point(CLUSTER,j)[i]-m[i];
			diff = diff*diff;
			v[i] += diff;
		}
		v[i] /= pts;

	}

	//store/return the variance
	//CLUSTER.VAR = v[i];
	//return v[i];
	CLUSTER.VAR = v;
	return v;
}

/********************************************************************[CLUSTER]*
 * Computes the CLUSTER's standard deviation (vector)
 * CLUSTER is the cluster that you want to know SD
 */
function cluster_calc_standard_deviation(CLUSTER){
	// The standard deviation is just the square-root of the variance
	var V = cluster_calc_variance(CLUSTER);
	var d = cluster_get_dimension(CLUSTER);
	var i; var SD = [];
	for(i=0; i<d; i+=1){
		SD[i] = Math.sqrt(V[i]);
	}
	CLUSTER.SD = SD;
	return SD;
}

/********************************************************************[CLUSTER]*
 * Returns a cluster's last computed center
 * CLUSTER is the cluster that you want to know the center of
 * call 'cluster_calc_center(CLUSTER)' if you want the updated center
 */
function cluster_get_center(CLUSTER){
	// DON'T recalculate the center here!
	// Certain algorithms benefit from controlling when the center is updated
	return CLUSTER.CENTROID;
}

/********************************************************************[CLUSTER]*
 * Merges CB into CA, ie appends points from B into A
 */
function cluster_merge(CA,CB){
	// Make sure we're merging clusters of the same dimensionallity
	if(CA.DIM !== CB.DIM){ return; }

	// Add all of the points of B to A
	var PTS = CB.POINTS;
	var i; var l = PTS.length;
	for(i=0; i<l; i+=1){
		cluster_add_point(CA,PTS[i]);
	}
}

/*
var T = cluster_create(2);
cluster_add_point(T,[0,0]);
cluster_add_point(T,[-1,-1]);
//console.info(cluster_get_point(T,0));
//console.info(cluster_get_point(T,1));
cluster_add_point(T,[-1,-1,0]);
//console.info(cluster_get_center(T));
var S = cluster_create(2);
cluster_add_point(S,[1,1]);
cluster_add_point(S,[1,-1]);
cluster_merge(S,T);
//console.info(cluster_get_center(S));
//console.info(S);
*/

/********************************************************************[K-MEANS]*
 * group_by_k_means(k,points) runs the k-means algorithm
 * It takes an array of points (a point is an array of numbers ie [0,1,...])
 * and returns an array of k CLUSTER objects that represent how those points
 * can be grouped
 */
function group_by_k_means(k,points){
	// We cannot have more clusters than points
	if(k > points.length){ k = points.length; }
	var dim = points[0].length;
	
	// Put the first k points into different clusters and compute the center
	var clusters = []; var centers = [];
	var i; var j;
	for(j=0; j<k; j+=1){
		clusters[j] = cluster_create(dim);
		cluster_add_point(clusters[j],points[j]);
		centers[j] = cluster_calc_center(clusters[j]);
	}
	
	// Limit the maximum number of iterations to 1000
	var maxIter = 1000;
	// If the mean hasn't changed we can stop the algorithm
	var isMeanNew = 1;

	var n = points.length;
	var min = 100000, minMean = 0, minIndex = 0, dv = [], dm;
	while(isMeanNew && --maxIter){
		
		// Clear the points stored in each group
		for(j=0; j<k; j+=1){
			cluster_clear_points(clusters[j]);
		}

		// Put each point into the 'closest' group
		n = points.length;
		for(i=0; i<n; i+=1){
			
			// Search for the closest group to this point
			min = 1000000; minMean = 0;
			for(j=0; j<k; j+=1){
				dv = vector_add(points[i],centers[j],-1);
				dm = vector_magnitude_squared(dv);
				if(dm < min){
					min      = dm;
					minIndex = j;
				}
			}

			cluster_add_point(clusters[minIndex],points[i]);
		}
		
		// Update the centers of each group
		// If any changed we'll need to iterate
		isMeanNew = 0;
		for(j=0; j<k; j+=1){
			dv = cluster_calc_center(clusters[j]);
			if(vector_is_equal(dv,centers[j])===0){
				isMeanNew = 1;
				centers[j] = dv;
			}
		}
	}

	return clusters;
}
//console.info(group_by_k_means(2,[[0,0],[1,1],[10,10],[-10,-10]]));

/***************************************************************[hierarchical]*
 * group_by_heirarchy(k,points,distFunc) heirarchically clusters the points
 * It takes an array of points (a point is an array of numbers (ie row vector))
 * and returns an array of k clusters that represent how those points can be grouped
 * distFunc is a function that computes the distance between the clusters
 * it should have the form distFunc(cluserA,clusterB)
 * where clusterA and clusterB are cluster objects
 */
function group_by_heirarchy(k,points,distance){
	// We cannot have more groups than points
	if(k > points.length){ k = points.length; }
	var dim = points[0].length;
	
	// Make each point it's own cluster
	var clusters = [];
	var i; var j; var l = points.length;
	for(i=0; i<l; i+=1){
		clusters[i] = cluster_create(dim);
		cluster_add_point(clusters[i],points[i]);
	}

	//n will be the count of clusters that we have
	var n = clusters.length;
	
	// We need to remember the minimum distance (and where)
	var d, minDis, im, jm, tempC;
	var maxIter = l;

	// While we have more clusters than we want
	while(n > k && maxIter--){

		// Start the minimum distance with data from our data set
		minDis = distance(clusters[0],clusters[1]);
		im = 0; jm = 1;

		// Find the minimum distance between between any two clusters
		// Note: we can start j at i+1 because:
		// (1) we don't want to compare the same cluster (i===j)
		// (2) i and j are indexing the same data so i=0 has handled
		// all of the comparisons for element 0 and one comparison for
		// elements 1..n; therefore, we don't need to compare anything
		// with element 0 again.
		for(i=0; i<n; i+=1){
			for(j=i+1; j<n; j+=1){

				//compute the distance
				d = distance(clusters[i],clusters[j]);

				// Remember the minimum
				if(d < minDis){
					minDis = d;
					im = i;
					jm = j;
				}
			}
		}

		// Merge the two closest clusters
		// Take out one of the clusters
		tempC = clusters.splice(jm,1)[0];

		//if the j cluster is before the i cluster in the array
		//then when the j cluster was removed everything after it
		//shifted down by 1 so update accordingly
		if(jm < im){im -= 1;}

		//put all of the elements of the removed custer into the other cluster
		cluster_merge(clusters[im],tempC);

		//clusters has now been updated so get the new count
		n = clusters.length;
	}

	//we've got k clusters (or took too long) so return the clusters
	return clusters;
}

/* This data set shows the difference between the distance methods for heir clusters
27,18
7,53
71,34
17,76
95,38
27,101
124,38
95,114
120,101
144,117
143,91
71,63
53,115
181,51
*/
/*******************************************************************[DISTANCE]*
 *
 */
function distance_nearest_neighbor(cA,cB){

	// Determine sizes
	var i; var Al = cluster_get_num_points(cA);
	var j; var Bl = cluster_get_num_points(cB);
	var dv; var dm;

	// Start the min with something in our data set
	var min = vector_magnitude_squared(vector_add(cluster_get_point(cA,0),cluster_get_point(cB,0),-1));

	// Find the minimum distance from every point in cA to every point in cB
	for(i=0; i<Al; i+=1){
		for(j=0; j<Bl; j+=1){
			dv = vector_add(cluster_get_point(cA,i),cluster_get_point(cB,j),-1);
			dm  = vector_magnitude_squared(dv); // don't waste sqrt
			if(dm < min){
				min = dm;
			}
		}
	}

	return min;
}
/*******************************************************************[DISTANCE]*
 *
 */
function distance_farthest_neighbor(cA,cB){

	// Determine sizes
	var i; var Al = cluster_get_num_points(cA);
	var j; var Bl = cluster_get_num_points(cB);
	var dv; var dm;

	// Start the max with something in our data set
	var max = vector_magnitude_squared(vector_add(cluster_get_point(cA,0),cluster_get_point(cB,0),-1));

	// Find the minimum distance from every point in cA to every point in cB
	for(i=0; i<Al; i+=1){
		for(j=0; j<Bl; j+=1){
			dv = vector_add(cluster_get_point(cA,i),cluster_get_point(cB,j),-1);
			dm  = vector_magnitude_squared(dv);
			if(dm > max){
				max = dm;
			}
		}
	}

	return max;
}
/*******************************************************************[DISTANCE]*
 *
 */
function distance_average_neighbor(cA,cB){

	// Determine sizes
	var i; var Al = cluster_get_num_points(cA);
	var j; var Bl = cluster_get_num_points(cB);
	var dv;

	//start the sum as 0
	var sum = 0;

	// Find the average distance from every point in cA to every point in cB
	for(i=0; i<Al; i+=1){
		for(j=0; j<Bl; j+=1){
			dv = vector_add(cluster_get_point(cA,i),cluster_get_point(cB,j),-1);
			sum += vector_magnitude_squared(dv);
		}
	}

	// Divide the sum by the number of times we added stuff to it
	// to get the average
	return sum/(Al*Bl);
}
/*******************************************************************[DISTANCE]*
 *
 */
function distance_center_neighbor(cA,cB){
	// Find the centers
	var midA = cluster_calc_center(cA);
	var midB = cluster_calc_center(cB);

	// Return the squared magnitude of the difference
	var dV = vector_add(midA,midB,-1);
	return vector_magnitude_squared(dV);
}

/*
//console.info(group_by_heirarchy(2,[[0,0],[1,1],[10,10],[-10,-10]],distance_center_neighbor));
//console.info(group_by_heirarchy(2,[[0,0],[1,1],[10,10],[-10,-10]],distance_nearest_neighbor));
//console.info(group_by_heirarchy(2,[[0,0],[1,1],[10,10],[-10,-10]],distance_farthest_neighbor));
//console.info(group_by_heirarchy(2,[[0,0],[1,1],[10,10],[-10,-10]],distance_average_neighbor));
*/
/********************************************************************[CLUSTER]*
 * This is a recursive clustering algorithm
 * cluster is the cluster that pt should be added to
 * points is an array of all points in the data set
 * d is the distance that allows points to be added to a cluster
 *  ie if the distance between any 2 points is <= d, then
 *     the points should be in the same cluster
 */
function group_by_clustering(cluster,pt,points,d){

	//if the point is already in a cluster, exit
	if(pt.C !== 0){return 0;}

	// Add the point to this cluster (and mark that it's in a cluster)
	cluster_add_point(cluster,pt.DATA);
	pt.C = 1;

	// Search for all the points that are close enough to this point
	var dV; var dM;
	var i=0; var l = points.length;
	for(i=0; i<l; i+=1){

		// Compute the distance 
		dV = vector_add(points[i].DATA,pt.DATA,-1);
		dM = vector_magnitude_squared(dV);

		// If the point is close enough, group it with this
		if(Math.sqrt(dM) <= d){
			// Add this point to the cluster.
			// This causes the newly added point to check 
			// if there are points close enough to include
			group_by_clustering(cluster,points[i],points,d);
		}
	}

	return 1;
}

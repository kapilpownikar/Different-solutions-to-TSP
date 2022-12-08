// Project Identifier: 9B734EC0C043C5A836EA0EBE4BEFEA164490B2C7

#include <iostream>
#include <iomanip>
#include <algorithm>
#include <stdio.h>
#include <string>
#include <vector>
#include <cmath>
#include <limits>

using namespace std;

struct obj {
	int x;
	int y;
	char type;
};

// Holds all the prim values for the given vertex
struct prim {
	char Kv = 'n'; // has vertex v been visited yet: 'v' = visited; 'n' = not visited
	double Dv = numeric_limits<double>::infinity(); // minimal edge weight to vertex V; root Dv set to 0
	int Pv = 0; // index of the vertex preceeding this vertex
};

struct arb {
	char Kv = 'n';
	double Dv = numeric_limits<double>::infinity();
};

// std euclidean distance func for part B and C
double euclidean(obj& a, obj& b) {

	int x1 = a.x; int x2 = b.x;
	int y1 = a.y; int y2 = b.y;

	return sqrt((double((x2 - x1)) * double((x2 - x1))) + (double((y2 - y1)) * double((y2 - y1))));
}

// non-square root one to save time for part C dist matrix
double euc(obj& a, obj& b) {
	int x1 = a.x; int y1 = a.y;
	int x2 = b.x; int y2 = b.y;

	return (double((x2 - x1)) * double((x2 - x1))) + (double((y2 - y1)) * double((y2 - y1)));
}

class OPT_TSP {
public:
	vector<int> path; // path
	int num_vertices; // total number of vertices in the graph
	vector<obj> vertices; // total num vertices in graph
	double best_dist = 0.00; // best distance so far
	vector<int> best_path; // best path we have seen so far
	double curr_dist_total = 0.00; // curr running total for our distance
	vector<vector<double>> matrix; // distance matrix


	// used separately for part C
	void read_vertices() {

		vertices.reserve(num_vertices);
		for (int i = 0; i < num_vertices; ++i) {
			int x; int y;
			cin >> x >> y;
			obj new_in = { x, y,'n' }; // obj used, char is irrelevant for this part (no need for visited unvisitied)
			vertices.push_back(new_in);
		}
	}

	// MST for part C
	// Used in Promising
	double min_spanning_tree_opttsp(int primlength, int permLength) {

		double totalcost = 0.00; // returning total cost here
		prim def = { 'n', numeric_limits<double>::infinity(), 0 };
		vector<prim> prims(primlength, def);

		// set the starting distance to 0
		prims[0].Dv = 0;

		for (int out = 0; out < int(path.size() - permLength); ++out) {

			// Step 1: Find the smallest Kv
			double shortest_dist = numeric_limits<double>::infinity();
			int Kv_idx = -1;
			for (int in = 0; in < int(path.size() - permLength); ++in) {

				// if it hasn't been discovered, find shortest of them all
				if (prims[in].Kv == 'n') {
					double dist = prims[in].Dv;
					if (dist < shortest_dist) {
						shortest_dist = dist;
						Kv_idx = in;
					}
				}
			} // end of loop we have our Kv and shortest distance

			// Step 2: Mark Kv as true; add to the running total
			prims[Kv_idx].Kv = 'v';
			totalcost += prims[Kv_idx].Dv;

			// Step 3: Loop over all vertices and update the false neighbors of Kv
			for (int i = 0; i < int(path.size() - permLength); ++i) {

				if (prims[i].Kv == 'n') {
					
					int idx_vertex = path[permLength + i];
					int idx_kv_vertex = path[permLength + Kv_idx];
					double distance = sqrt(matrix[idx_vertex][idx_kv_vertex]);
					// for each vertex w adjacent to v, check whether Dw is greater than curr_dist,
					// if it is, then new Dw = curr_dist, and new Pw = v
					if (prims[i].Dv > distance) {
						prims[i].Dv = distance;
						prims[i].Pv = Kv_idx;
					}
				}

			} // end Step 3

		} // End of looping through all vertices

		return totalcost;
	}

	bool promising(int permLength) {

		bool promise = false;
		//// O(n!) is < O(n^2 + 2n) for size diff less than 4 so we might as well just explore those branches
		if (int(path.size() - permLength) < 5) {
			return true;
		}

		double mst_cost = min_spanning_tree_opttsp(int(path.size() - permLength), permLength);

		double arm1 = numeric_limits<double>::infinity();
		double arm2 = numeric_limits<double>::infinity();
		
		for (int i = 0; i < int(path.size() - permLength); ++i) {
			
			double curr_arm1 = sqrt(matrix[path[permLength + i]][path[0]]); // start to potential
			double curr_arm2 = sqrt(matrix[path[permLength + i]][path[permLength - 1]]); // end to potential
			if (curr_arm1 < arm1) {
				arm1 = curr_arm1;
			}
			if (curr_arm2 < arm2) {
				arm2 = curr_arm2;
			}
		} // at the end we'll have the min distance to connect both arms

		double potential_cost = mst_cost + arm1 + arm2;

		double wt_projected = potential_cost + curr_dist_total; // partial soln cost + projected cost with mst

		if (wt_projected > best_dist) {
			promise = false;
		}
		else {
			promise = true;
		}

		//for (size_t i = 0; i < path.size(); ++i)
		//	cout << setw(2) << path[i] << ' ';
		//cout << setw(4) << permLength << setw(10) << curr_dist_total;
		//cout << setw(10) << arm1 << setw(10) << arm2;
		//cout << setw(10) <<  mst_cost << setw(10) << wt_projected << "  " << boolalpha << promise << '\n';

		return promise;
	}

	void genPerms(uint32_t permLength) {
		if (permLength == path.size()) {
			// Do something with the path
			double edge_last = sqrt(matrix[0][path[path.size()-1]]);
			curr_dist_total += edge_last;
			if (curr_dist_total < best_dist) {
				best_dist = curr_dist_total;
				best_path = path;
				/*cout << "New best cost achieved: " << best_dist << "\n";*/
			}
			curr_dist_total -= edge_last;
			return;
		}  // if ..complete path

		if (!promising(permLength)) {
			return;
		}  // if ..not promising

		for (size_t i = permLength; i < path.size(); ++i) {
			swap(path[permLength], path[i]);
			double curr_edge = sqrt( matrix[path[permLength-1]][path[permLength]] );
			curr_dist_total += curr_edge;
			genPerms(permLength + 1);
			curr_dist_total -= curr_edge;
			swap(path[permLength], path[i]);
		}  // for ..unpermuted elements
	}  // genPerms()

	double upper_bound_tsp(vector<vector<double>> &matrix) {

		prim def_ault = { 'n', numeric_limits<double>::infinity(), 0 };
		vector<prim> fastprim(num_vertices, def_ault);
		double total_tour = 0.00;

		// variables used in looping
		double shortest_itok = 0.00; double shortest_jtok = 0.00;
		int vertex_curr = 0;

		// Step 1: start by adding 2 vertices to the path (0 and 1 idx)
		path.push_back(0); path.push_back(0);
		// fastprim[0].Kv = 'v'; fastprim[1].Kv = 'v';
		fastprim[0].Dv = 0; fastprim[1].Dv = 0;
		fastprim[0].Pv = 0; fastprim[1].Pv = 0;

		total_tour += fastprim[0].Dv + fastprim[1].Dv; // Adding this to the total tour

		for (int outer = 1; outer < num_vertices; ++outer) {

			double shortest_so_far = numeric_limits<double>::infinity(); // shortest distance so far starts as inf

			// loops through all nodes in the path 
			// path circle should also consider last and first as i and j
			// Eg: path: 0 -> 1 -> 2 -> 3 ; Here: pairs of i & j are: { (0,1), (1,2), (2,3), (3,0) }
			for (size_t vtx = 0; vtx < path.size() - 1; ++vtx) {

				double dist_ij = sqrt(matrix[path[vtx]][path[vtx + 1]]);
				double dist_itok = sqrt(matrix[path[vtx]][outer]);
				double dist_jtok = sqrt(matrix[path[vtx + 1]][outer]);

				double diff_eval = (dist_itok + dist_jtok) - dist_ij; // this is the C(i,k) + C(j,k) - C(i,j)

				if (diff_eval < shortest_so_far) {

					shortest_so_far = diff_eval;
					shortest_itok = dist_itok;
					shortest_jtok = dist_jtok;
					vertex_curr = int(vtx);
				}

			} // End inner loop

			path.insert((path.begin() + vertex_curr + 1), outer); // insert 'arbitrary' new vertex chosen
			total_tour += shortest_so_far; // add to total length

			// update values for new added element
			fastprim[outer].Pv = path[vertex_curr]; // check if % needed
			fastprim[outer].Dv = shortest_itok; // this is all w.r.t i as prev since i -> k -> j when inserted
			fastprim[outer].Kv = 'v';

			fastprim[ path[int(vertex_curr + 2) % int(path.size())] ].Dv = shortest_jtok;
			fastprim[ path[int(vertex_curr + 2) % int(path.size())] ].Pv = outer;
		} // End outer loop
		path.pop_back();

		return total_tour;
	}

	// Optimal TSP configuration
	// Uses Branch-and-Bound method
	// NOT EFFICIENT in time and memory to >40 vertices since it uses a distance matrix
	void opt_tsp() {

		read_vertices(); // reserves vertices capacity and reads coordinates

		// populating distance matrix first
		matrix.reserve(num_vertices);
		for (size_t row = 0; row < vertices.size(); ++row) {
			vector<double> vect(num_vertices, 0.00);
			for (size_t col = 0; col < vertices.size(); ++col) {
				double dist = euc(vertices[row], vertices[col]);
				vect[col] = dist;
			}
			matrix.push_back(vect);
		}

		// Step 1: Calculate the upperbound and set it to best dist so far; also changes path vector to fasttsp path
		// Uses FASTTSP heuristic for this part
		best_dist = upper_bound_tsp(matrix);
		best_path = path;

		// resetting path for debug output to 0->1->2->3....
		/*for (size_t i = 0; i < path.size(); ++i) {
			path[i] = int(i);
		}*/

		// Step 2: Call genperms on 1
		genPerms(1);
		// ----- // ----- // ----- // ----- //
		cout << best_dist << "\n";

		for (size_t i = 0; i < path.size(); ++i) {
			cout << best_path[i] << " ";
		}
		return;
	}
};


class GameADT {
public:

	vector<obj> vertices = {};
	int num_vertices;
	vector<int> path = {}; // used for both fasttsp and opttsp

	char modetype; // [MST] {FASTTSP] [OPTTSP} 
	obj starting = {0, 0, 'D'}; // starting point at (0,0)

	// Vertex creation function
	obj create_obj(int& x, int& y) {

		if ( (x == 0 && y < 0) 
			|| (x < 0 && y == 0) || (x == 0 && y == 0) ) {
			return obj{ x, y, 'D' }; // vertex is on decontamination zone
		}

		if (x < 0 && y < 0) {
			return obj{ x, y, 'L' }; // vertex is in Lab
		}

		return obj{ x, y, 'O' }; // vertex is in the outer regions
	}

	// MST Distance: Euclidean Distance
	// Assumption: Dist b/w Outer regions or Outer to Decontam is Euclidean
	// Outer to Inner or Inner to Outer == Infinity 
	double mst_dist(obj& pos_a, obj& pos_b) {

		// local dist variables
		int x1 = pos_a.x; int x2 = pos_b.x;
		int y1 = pos_a.y; int y2 = pos_b.y;
		double base = 0.0;
		char atype = pos_a.type; char btype = pos_b.type;

		if (atype == 'L' && btype == 'O') {
			return numeric_limits<double>::infinity();
		}
		else if (btype == 'L' && atype == 'O') {
			return numeric_limits<double>::infinity();
		}
		else {
			// Normal Euclidean Distance return 
			// TODO: Check if this works
			return sqrt((double((x2 - x1)) * double((x2 - x1))) + (double((y2 - y1)) * double((y2 - y1))));
		}

		return base; // default value to handle control of non-void function error
	}

	// part 1: MST Implementation
	// - Objects needed: - Euclidean Distance
	//					 - Prim's values
	void Min_Spanning_Tree() {
		
		vector<prim> primvals = {}; // set in only when MST mode required

		vertices.reserve(num_vertices); // reserving capacity
		primvals.reserve(num_vertices); // reserving capacity
		bool out_zone = false; bool lab_zone = false; bool decontam = false; // used to check error condition

		for (int i = 0; i < num_vertices; ++i) {
			int x_in; int y_in;
			cin >> x_in >> y_in;
			obj ob_in = create_obj(x_in, y_in);
			prim primobj = { 'n', numeric_limits<double>::infinity(), 0 };
			// Triggers only the 1st time
			if (!lab_zone && ob_in.type == 'L') { lab_zone = true; } 
			else if (!decontam && ob_in.type == 'D') { decontam = true; }
			else if (!out_zone && ob_in.type == 'O') { out_zone = true; }

			vertices.push_back(ob_in); // adds vertices into the vector
			primvals.push_back(primobj); // adds default prims into the vector 
		}

		// check first for error conditon
		if (num_vertices > 2) {
			if (!decontam && (lab_zone && out_zone)) {
				cerr << "Cannot construct MST\n";
				exit(1);
			}
		}

		double total_dist = 0.00; // this value increments
		// Starting point is first vertex
		primvals[0].Dv = 0;

		for (int out = 0; out < num_vertices; ++out) {

			// Step 1: find the smallest Kv
			double shortest_dist = numeric_limits<double>::infinity();
			int closest = -1;
			for (int in = 0; in < num_vertices; ++in) {

				// if it hasn't been discovered
				if (primvals[in].Kv == 'n') {
					double dist = primvals[in].Dv;
					if (dist < shortest_dist) {
						shortest_dist = dist;
						closest = in;
					}
				}
			} // At the end of this inner loop we have a shortest index and a shortest distance

			// Step 2: Mark Kv as true; add to running total
			primvals[closest].Kv = 'v';
			total_dist += primvals[closest].Dv;

			// Step 3: Loop over all vertices and update false neighbours of Kv
			for (int i = 0; i < num_vertices; ++i) {

				if (primvals[i].Kv == 'n') {

					double curr_dist = mst_dist(vertices[closest], vertices[i]);
					// for each vertex w adjacent to v, check whether Dw is greater than curr_dist, 
					// if it is, then new Dw = curr_dist, and new Pw = v
					if (primvals[i].Dv > curr_dist) {
						primvals[i].Dv = curr_dist;
						primvals[i].Pv = closest;
					}

				}
			} // End of step 3

		} // end of looping through all vertices

		cout << total_dist << "\n";

		// Starts with index 1 since idx 0 was the root
		for (int j = 1; j < num_vertices; ++j) {

			// Only if visited
			if (primvals[j].Kv == 'v') {
				int pred_j = primvals[j].Pv;
				if (j < pred_j) { cout << j << " " << pred_j << "\n"; }
				else { cout << pred_j << " " << j << "\n"; }
			}
		}

	}

	void print_fast_tsp(double & total_tour, vector<int> &path) {

		// Printing the right output
		cout << total_tour << "\n";
		// loop to print path
		for (size_t i = 0; i < path.size(); ++i) {
			cout << path[i] << " ";
		}
		cout << "\n";

	}

	// --mode FASTTSP
	// Uses Arbitrary Insertion Heuristic
	// Disregards zones so all objects belong to 'O' (makes it easier)
	void fast_tsp() {

		path.reserve(num_vertices + 1);
		vector<prim> fastprim; fastprim.reserve(num_vertices + 1);
		vertices.reserve(num_vertices);
		double total_tour = 0.00;

		for (int i = 0; i < num_vertices; ++i) {
			int x_in, y_in;
			cin >> x_in >> y_in;
			obj obnew = { x_in, y_in, 'n' };
			vertices.push_back(obnew);
			prim prnew = { 'n', numeric_limits<double>::infinity(), -1 }; // -1 could be a cause for vector subscript out of range
			fastprim.push_back(prnew);
		}

		// variables used in the loops
		double shortest_itok = 0.00; double shortest_jtok = 0.00;
		int vertex_curr = 0;

		// Step 1: start by adding 2 vertices to the path (0 and 1 idx)
		path.push_back(0); path.push_back(0);
		// fastprim[0].Kv = 'v'; fastprim[1].Kv = 'v';
		fastprim[0].Dv = 0; fastprim[1].Dv = 0;
		fastprim[0].Pv = 0; fastprim[1].Pv = 0;

		total_tour += fastprim[0].Dv + fastprim[1].Dv; // Adding this to the total tour

		// Step 2: looping condition here:
		// the outer variable acts as the arbitrary node to insert;
		// not actually arbitrarily selected
		for (int outer = 1; outer < num_vertices; ++outer) {

			double shortest_so_far = numeric_limits<double>::infinity(); // shortest distance so far starts as inf

			// loops through all nodes in the path 
			// path circle should also consider last and first as i and j
			// Eg: path: 0 -> 1 -> 2 -> 3 ; Here: pairs of i & j are: { (0,1), (1,2), (2,3), (3,0) }
			for (size_t vtx = 0; vtx < path.size() - 1; ++vtx) {

				double dist_ij = euclidean(vertices[path[vtx]], vertices[path[vtx + 1]]);
				double dist_itok = euclidean(vertices[path[vtx]], vertices[outer]);
				double dist_jtok = euclidean(vertices[path[vtx + 1]], vertices[outer]);

				double diff_eval = (dist_itok + dist_jtok) - dist_ij; // this is the C(i,k) + C(j,k) - C(i,j)

				if (diff_eval < shortest_so_far) {

					shortest_so_far = diff_eval;
					shortest_itok = dist_itok;
					shortest_jtok = dist_jtok;
					vertex_curr = int(vtx);
				}
				
			} // End inner loop

			path.insert((path.begin() + vertex_curr + 1), outer); // insert 'arbitrary' new vertex chosen
			total_tour += shortest_so_far; // add to total length

			// update values for new added element
			fastprim[outer].Pv = path[vertex_curr]; // check if % needed
			fastprim[outer].Dv = shortest_itok; // this is all w.r.t i as prev since i -> k -> j when inserted
			// fastprim[outer].Kv = 'v';

			fastprim[path[(vertex_curr + 2) % path.size()]].Dv = shortest_jtok;
			fastprim[path[(vertex_curr + 2) % path.size()]].Pv = outer;
		} // End outer loop
		
		path.pop_back(); // remove the final zero added

		if (modetype == 'F') {
			print_fast_tsp(total_tour, path);
		}
	}

};
#include <algorithm>
#include <climits>
#include <fstream>
#include <iostream>
#include <map>
#include <mutex>
#include <random>
#include <set>
#include <stdexcept>
#include <sstream>
#include <thread>
#include <utility>
#include <vector>

using std::cin;
using std::cout;
using std::endl;
using std::ifstream;
using std::istream;
using std::logic_error;
using std::make_pair;
using std::map;
using std::max;
using std::min;
using std::mt19937;
using std::mutex;
using std::ofstream;
using std::ostream;
using std::pair;
using std::random_device;
using std::random_shuffle;
using std::runtime_error;
using std::sort;
using std::stable_partition;
using std::string;
using std::stringstream;
using std::set;
using std::thread;
using std::uniform_int_distribution;
using std::vector;

typedef unsigned int vertex_id;

// Constants
static int NUM_SIMULATIONS = 300;
static const int NUM_THREADS = 8;
static const bool SHUFFLE = false;
static const bool DEBUG_MESSAGES = false;

// Locks
static mutex ostream_lock;
static mutex statistics_lock;

// Random number generation variables
static random_device rd;
static mt19937 mt( rd() );

// ============================================================================
// ============================ CLASS DEFINITIONS =============================
// ============================================================================

class Simple_Graph {
public:
	const map<vertex_id, set<vertex_id> >& edges() const { return graph; }
    
	void clear() { graph.clear(); }

    void clear_edges();
    
    void add_vertex(vertex_id id, unsigned int weight);

    void add_edge(vertex_id a, vertex_id b);

	void set_weight(vertex_id id, unsigned int weight);

	unsigned int optimal_vertex_cover(set<vertex_id>& cover) const;

	unsigned int greedy_vertex_cover(set<vertex_id>& cover) const;

	unsigned int pricing_vertex_cover(set<vertex_id>& cover) const;

	unsigned int improved_pricing_vertex_cover(set<vertex_id>& cover) const;

    string to_string() const;

	string to_string_output() const;

    void generate_as_random(unsigned int num_vertices, unsigned int num_edges,
        unsigned int max_weight);


private:
    map<vertex_id, set<vertex_id> > graph;
    map<vertex_id, unsigned int> weights;

	bool is_valid_cover(const set<vertex_id> &combination) const;

	bool equal_weights(const pair<vertex_id, vertex_id>& edge) const;

    vertex_id greatest_ratio_vertex(
        const map<vertex_id, set<vertex_id> > &graph) const;

	unsigned int pricing_vertex_cover( set<vertex_id>& cover, bool improved ) const;

    vertex_id pick_random(const set<vertex_id>& vertex_set) const;
    
    void add_random_edges(const set<vertex_id>& vertex_set, 
        unsigned int num_new_edges);
    
    vector<vertex_id> loop_erasing_random_walk(
        const set<vertex_id>& tree_vertices, const set<vertex_id>& vertex_set) const;

    void wilson_spanning_tree(const set<vertex_id>& vertex_set);
};

// ============================================================================



// ============================================================================
// ====================== MISCELLANEOUS GRAPH FUNCTIONS =======================
// ============================================================================

string Simple_Graph::to_string() const
{
    stringstream ss;
    for (map<vertex_id, set<vertex_id> >::const_iterator i = graph.begin();
        i != graph.end(); i++)
    {
        vertex_id current = i->first;
        unsigned int weight = weights.at(current);
        const set<vertex_id>& adjacent = i->second;

        ss << current << " (w. " << weight << ")" << ": ";

        for(set<vertex_id>::const_iterator j = adjacent.begin();
            j != adjacent.end(); j++) 
        {
            if ( j != adjacent.begin() ) 
				ss << ", ";
            ss << *j;
        }
		ss << endl;
    }

    return ss.str();
}

string Simple_Graph::to_string_output() const
{
	stringstream ss;

	ss << graph.size() << endl;

	for (map<vertex_id, set<vertex_id> >::const_iterator i = graph.begin();
		i != graph.end(); i++)
	{
		vertex_id current = i->first;
		const set<vertex_id>& adjacent = i->second;

		ss << weights.at(i->first);

		for (set<vertex_id>::const_iterator j = adjacent.begin();
			j != adjacent.end(); j++)
		{
			ss << " " << *j;
		}
		ss << " 0" << endl;
	}

	return ss.str();
}

vertex_id Simple_Graph::greatest_ratio_vertex(
    const map<vertex_id, set<vertex_id> > &graph) const
{
    double best_ratio = -1;
    vertex_id best_vertex;
    typedef map<vertex_id, set<vertex_id> >::const_iterator vertex_iterator;
    
	if(DEBUG_MESSAGES)
		cout << "Vertices:" << endl;

    for(vertex_iterator v_it = graph.begin(); v_it != graph.end(); v_it++)
    {
        vertex_id current_vertex = v_it->first;
        set<vertex_id>::size_type degree = v_it->second.size();

		double vertex_weight = weights.find(current_vertex)->second;
		double ratio = degree / vertex_weight;

		if (DEBUG_MESSAGES)
		{
			cout << "    Vertex: " << current_vertex << endl;
			cout << "    Degree: " << degree << endl;
			cout << "    Weight: " << vertex_weight << endl;
			cout << "    Ratio: " << ratio << endl;
			cout << endl;
		}

        if ( ratio > best_ratio )
        {
            best_ratio = ratio;
            best_vertex = current_vertex;
        }
    }

    return best_vertex;
}

bool Simple_Graph::equal_weights(const pair<vertex_id, vertex_id>& edge) const
{
    vertex_id vertex1 = edge.first;
    vertex_id vertex2 = edge.second;

    return weights.at(vertex1) == weights.at(vertex2);
}

bool Simple_Graph::is_valid_cover(const set<vertex_id> &combination) const
{
	for (map<vertex_id, set<vertex_id> >::const_iterator it = graph.begin();
		it != graph.end(); it++)
	{
		vertex_id current_vertex = it->first;
		const set<vertex_id> &edges = it->second;

		// Continue if the current vertex is in combination
		if (combination.find(current_vertex) != combination.end())
			continue;

		for (set<vertex_id>::const_iterator edge_it = edges.begin();
			edge_it != edges.end(); edge_it++)
		{
			vertex_id other_vertex = *edge_it;

			// Check if the second vertex is in combination
			if (combination.find(other_vertex) == combination.end())
				return false;   // Found an edge not covered           
		}
	}

	return true;
}


// ============================================================================



// ============================================================================
// ========================= GRAPH BUILDING FUNCTIONS =========================
// ============================================================================

void Simple_Graph::clear_edges()
{
    for(map<vertex_id, set<vertex_id> >::iterator it = graph.begin();
        it != graph.end(); it++)
    {
        it->second.clear();
    } 
}

void Simple_Graph::add_vertex(vertex_id id, unsigned int weight)
{
    // Check if the vertex has already been added
    if ( graph.find(id) != graph.end() )
    {
        stringstream ss;
        ss << "Vertex " << id << " already in the graph" << endl;
        throw logic_error( ss.str() );
    }

    if ( weight <= 0 )
    {
        stringstream ss;
        ss << "Vertex " << id << " must have a positive weight" << endl;
        throw logic_error( ss.str() );
    }

    // Create the vertex in the graph, with no edges
    graph[id];
    weights[id] = weight;
}

void Simple_Graph::add_edge(vertex_id a, vertex_id b)
{
    stringstream ss;

    // ===== Check if both vertices exist =====
    if ( graph.find(a) == graph.end() )
    {
        ss << "Invalid edge, vertex " << a << " does not exist.";
        throw logic_error( ss.str() );
    }
    
    if ( graph.find(b) == graph.end() )
    {
        ss << "Invalid edge, vertex " << b << " does not exist.";
        throw logic_error( ss.str() );
    }
    // ========================================

    // Insert the other vertex id in each vertex's adjacency set
    graph[a].insert(b);
    graph[b].insert(a);
}

void Simple_Graph::set_weight(vertex_id id, unsigned int weight)
{
	if (weights.find(id) == weights.end())
		throw logic_error("Cannot assign a weight to a non-existent vertex");

	weights[id] = weight;
}

vertex_id Simple_Graph::pick_random(const set<vertex_id>& vertex_set) const
{
    std::uniform_int_distribution<set<vertex_id>::size_type> dist( 0, vertex_set.size() - 1 );
    const set<vertex_id>::size_type rand_idx = dist(mt);

    // Advance it to the element at position rand_idx
    set<vertex_id>::const_iterator it = vertex_set.begin();
    for ( set<vertex_id>::size_type i = 0; i < rand_idx; i++ )
        it++;

    return *it;
}

vector<vertex_id> Simple_Graph::loop_erasing_random_walk(
    const set<vertex_id>& tree_vertices, const set<vertex_id>& vertex_set) const
{
    vector<vertex_id> walk;
    map<vertex_id, vertex_id> parent_map;

    // Check if there are vertices which are not part of the tree
    if ( tree_vertices.size() == vertex_set.size() )
        return walk;

    // Pick a random unused vertex as the start vertex
    vertex_id current = pick_random(vertex_set);

    while( tree_vertices.find(current) != tree_vertices.end() )
        current = pick_random(vertex_set);
    
    // Iterate until we have arrived at a vertex that is part of the tree
    while ( tree_vertices.find(current) == tree_vertices.end() )
    {
        // Pick another vertex to be the next in the walk
        vertex_id neighbour = pick_random(vertex_set);
        
        while( neighbour == current)
            neighbour = pick_random(vertex_set);

        /*  Check if the selected vertex is already in the walk. If so, we have
            found a loop and must remove it */
        if ( parent_map.find(neighbour) != parent_map.end() )
        {
            vertex_id last_valid = neighbour;

            // Go back through the walk, removing the loop
            while ( current != last_valid )
            {
                vertex_id parent = parent_map[current];
                parent_map.erase(current);
                current = parent;
            }
        }
        // The selected vertex is a new vertex, so we need only update the state
        else
        {
            parent_map[neighbour] = current;
            current = neighbour;
        }
    }

    // ===== Transfer the random walk to the return vector =====
    walk.resize( parent_map.size() + 1 );

    for (vector<vertex_id>::iterator rev_it = walk.begin();
        rev_it != walk.end(); rev_it++)
    {
        *rev_it = current;
        current = parent_map[current];
    }
    // =========================================================

    return walk;
}

void Simple_Graph::wilson_spanning_tree(const set<vertex_id>& vertex_set)
{
    // Create the set of vertices used in the spanning tree
    set<vertex_id> tree_vertices;
    // Add vertex 1 as the root
    tree_vertices.insert(1);

    vector<vertex_id> walk = loop_erasing_random_walk(tree_vertices, vertex_set);
    while( !walk.empty() )
    {
        tree_vertices.insert( walk.begin(), walk.end() );

        for (vector<vertex_id>::size_type i = 0; i < walk.size() - 1; i++)
            add_edge( walk[i], walk[i+1] );

        walk = loop_erasing_random_walk(tree_vertices, vertex_set);
    }
}

void Simple_Graph::add_random_edges(const set<vertex_id>& vertex_set, 
    unsigned int num_new_edges)
{
    /*  Store all the vertices that can still support new edges. A vertex will be 
        removed from this set if it already has |V|-1 edges, i.e., it has an edge 
        to every other vertex */
    set<vertex_id> available_vertices = vertex_set;

    while( num_new_edges > 0  && available_vertices.size() > 0)
    {
        // ===== Pick two different vertices =====
        vertex_id vertex1 = pick_random(available_vertices);
        
        /*  Compute the set of available vertices not adjacent to vertex1. 
            These are the options that can have an edge to vertex1 */
        set<vertex_id> not_adjacent = available_vertices;
        not_adjacent.erase( vertex1 );
        set<vertex_id>& vertex1_adjacent = graph[vertex1];
        
        for (set<vertex_id>::const_iterator it = vertex1_adjacent.begin();
            it != vertex1_adjacent.end(); it++)
        {
            // Check if the vertex *it is an available one
            set<vertex_id>::iterator pos = not_adjacent.find( *it );
            /*  If the vertex *it is no longer available, remove it from 
                the not_adjacent list */
            if ( pos != not_adjacent.end() )
                not_adjacent.erase( pos );
        }

        vertex_id vertex2 = pick_random(not_adjacent);

        // =======================================

        set<vertex_id>& vertex2_adjacent = graph[vertex2];
        
        vertex1_adjacent.insert(vertex2);
        vertex2_adjacent.insert(vertex1);

        if ( vertex1_adjacent.size() == vertex_set.size() - 1 )
            available_vertices.erase( vertex1 );
        if ( vertex2_adjacent.size() == vertex_set.size() - 1 )
            available_vertices.erase(vertex2);

        
        num_new_edges--;
    }
}

void Simple_Graph::generate_as_random(unsigned int num_vertices, 
    unsigned int num_edges, unsigned int max_weight)
{
    if (num_edges < num_vertices - 1)
    {
        stringstream ss;
        ss << "A graph with |V| = " << num_vertices << " and |E| = " << 
            num_edges << "cannot be connected" << endl;
        throw logic_error( ss.str() );
    }
 
    if (num_edges > num_vertices * (num_vertices - 1) / 2)
    {
        stringstream ss;
        ss << "A graph with |V| = " << num_vertices << " can have to most " << 
            (num_vertices * (num_vertices - 1) / 2) << " edges" << endl;
        throw logic_error( ss.str() );
    }

    clear();

    set<vertex_id> vertex_set;
    uniform_int_distribution<unsigned int> weight_dist(1, max_weight);

    // Add the specified number of vertices to both the graph and the vertex set
    for (unsigned int i = 1; i <= num_vertices; i++)
    {
        add_vertex(i, weight_dist(mt));
        vertex_set.insert(i);
    }

    // Generate a random spanning tree to ensure the graph is connected
    wilson_spanning_tree(vertex_set);

    /*  Add random edges to achieve the specified density. So far, there are
        num_vertices-1 edges, since the graph is a tree */
    add_random_edges(vertex_set, num_edges - (num_vertices-1));
}

// ============================================================================



// ============================================================================
// ========================= VERTEX COVER FUNCTIONS =========================
// ============================================================================

unsigned int Simple_Graph::optimal_vertex_cover(set<vertex_id> &cover) const
{
	cover.clear();
	unsigned int best_weight_sum = UINT_MAX;

	// This set will hold the vertex sets checked as possible covers
	set<set<vertex_id> > valid_covers;
	set<set<vertex_id> >* subsets = new set<set<vertex_id> >;
	set<set<vertex_id> >* next_subsets = new set<set<vertex_id> >;

	// Insert all sets of size 1 in subsets
	for (map<vertex_id, set<vertex_id> >::const_iterator it = graph.begin();
		it != graph.end(); it++)
	{
		set<vertex_id> unit_set;
		unit_set.insert(it->first);

		// If the unit set is a valid cover, return it
		if ( is_valid_cover(unit_set) )
		{
			cover = unit_set;
			delete subsets;
			delete next_subsets;
			return weights.at( it->first );
		}

		subsets->insert(unit_set);
	}

	while (true)
	{
		// Remember the min weight from all subsets in this iteration
		unsigned int min_iteration_weight = UINT_MAX;
		
		for (set<set<vertex_id> >::const_iterator comb_it = subsets->begin();
			comb_it != subsets->end(); comb_it++)
		{
			const set<vertex_id> &current_subset = *comb_it;


			/*  For each vertex and each subset, add a subset that contains the
			vertex and one that does not */
			for (map<vertex_id, set<vertex_id> >::const_iterator vertex_it = graph.begin();
				vertex_it != graph.end(); vertex_it++)
			{
				vertex_id current_vertex = vertex_it->first;

				// Continue if the current vertex is already in it
				if (current_subset.find(current_vertex) != current_subset.end())
					continue;

				set<vertex_id> current_subset_copy = current_subset;

				/*  Add the current subset U {current_vertex} to the next
				iteration set */
				current_subset_copy.insert(current_vertex);

				// Compute the weight sum of this solution
				unsigned int weight_sum = 0;
				for (set<vertex_id>::const_iterator cover_it = current_subset_copy.begin();
					cover_it != current_subset_copy.end(); cover_it++)
				{
					weight_sum += weights.at(*cover_it);
				}

				// Update the min subset weight in this iteration
				if (weight_sum < min_iteration_weight)
					min_iteration_weight = weight_sum;

				// Check if the subset is a better solution
				if (weight_sum < best_weight_sum && is_valid_cover(current_subset_copy))
				{
					cover = current_subset_copy;
					best_weight_sum = weight_sum;
				}

				next_subsets->insert(current_subset_copy);
			}
		}

		/*	Stop if no subsets in this iteraiton have a smaller weight sum than the 
			current solution. Since the subsets always grow from one iteration to
			the next, no future subsets are going to be better solutions. */
		if ( !cover.empty() && best_weight_sum <= min_iteration_weight )
			return best_weight_sum;

		// Update the subsets pointer
		delete subsets;
		subsets = next_subsets;
		next_subsets = new set<set<vertex_id> >;
	}

	delete subsets;
	delete next_subsets;
}

unsigned int Simple_Graph::greedy_vertex_cover(set<vertex_id>& cover) const
{
    cover.clear();

	// Counter for the total weight sum of the vertex cover
	unsigned int cover_weight_sum = 0;

    map<vertex_id, set<vertex_id> > graph_copy = graph;

    while( !graph_copy.empty() )
    {
        // Find the vertex with the greatest degree/weight ratio
        vertex_id grv = greatest_ratio_vertex(graph_copy);

		if (DEBUG_MESSAGES)
		{
			cout << "Greatest ratio vertex: " << grv << endl;
			cout << "********************" << endl << endl;
		}

		// Add the vertex to the vertex cover
        cover.insert( grv );
		cover_weight_sum += weights.at(grv);

        /*  Remove the greatest ratio vertex from the adjacency list of each
            of its adjacent vertices */
        set<vertex_id> &adj_vertices = graph_copy[grv];
        for (set<vertex_id>::const_iterator adj_vertex_it = adj_vertices.begin();
            adj_vertex_it != adj_vertices.end(); adj_vertex_it++)
        {
            vertex_id adj_vertex = *adj_vertex_it;
            set<vertex_id> &adj_vertex_edges = graph_copy[adj_vertex];
            
            adj_vertex_edges.erase( grv );
            
            // If the vertex has become isolated, remove it from graph_copy
            if ( adj_vertex_edges.empty() )
                graph_copy.erase( *adj_vertex_it );
        }

        // Remove the greatest degree vertex from the graph
        graph_copy.erase( grv );
    }

	return cover_weight_sum;
}

unsigned int Simple_Graph::pricing_vertex_cover(set<vertex_id>& cover) const
{
    return pricing_vertex_cover(cover, false); 
}

unsigned int Simple_Graph::improved_pricing_vertex_cover(set<vertex_id>& cover) const
{
    return pricing_vertex_cover(cover, true); 
}

unsigned int rand_uint(unsigned int range_limit)
{
	uniform_int_distribution<unsigned int> dist(0, range_limit - 1);
	return dist(mt);
}

unsigned int Simple_Graph::pricing_vertex_cover( set<vertex_id>& cover, bool improved) const
{
    typedef pair<vertex_id, vertex_id> edge;

    /*  Temporary set containing all edges. This is created to avoid duplicated 
        edges in the following edge vector */
    set<edge> edge_set;

    for (map<vertex_id, set<vertex_id> >::const_iterator graph_it = graph.begin();
        graph_it != graph.end(); graph_it++)
    {
        vertex_id current_vertex = graph_it->first;
        const set<vertex_id>& adjacents = graph_it->second;
        
        for (set<vertex_id>::const_iterator set_it = adjacents.begin();
            set_it != adjacents.end(); set_it++)
        {
            // Identify the incident vertices of smaller and greater id
            vertex_id adjacent_vertex = *set_it;
            vertex_id min_vertex = min(current_vertex, adjacent_vertex);
			vertex_id max_vertex = max(current_vertex, adjacent_vertex);
            
            // Add an edge to the edge set with the incident vertices ordered by id
            edge_set.insert( edge(min_vertex, max_vertex) );
        }
    }

    // Vector containing all the edges
    vector<edge> edge_vec( edge_set.begin(), edge_set.end() );

	if (DEBUG_MESSAGES)
	{
		cout << "Edge vector:" << endl;
		for (auto it = edge_vec.begin(); it != edge_vec.end(); it++)
			cout << "  (" << it->first << ", " << it->second << ")" << endl;
	}

	if (SHUFFLE)
	{
		random_shuffle(edge_vec.begin(), edge_vec.end(), rand_uint);

		if (DEBUG_MESSAGES)
		{
			cout << "Shuffled edge vector:" << endl;
			for (auto it = edge_vec.begin(); it != edge_vec.end(); it++)
				cout << "  (" << it->first << ", " << it->second << ")" << endl;
		}
	}

    /*  Improved pricing heuristic. This heuristic partitions the vector in such
        a way that all edges whose vertices have different weights will be chosen
        before edges whose vertices have equal weights */
    if (improved)
    {
        stable_partition(edge_vec.begin(), edge_vec.end(), 
            [this] (const pair<vertex_id, vertex_id>& edge) { return this->equal_weights(edge); } );
        
		if (DEBUG_MESSAGES)
		{
			cout << "Edge vector after heuristic:" << endl;
			for (auto it = edge_vec.begin(); it != edge_vec.end(); it++)
				cout << "  (" << it->first << ", " << it->second << ")" << endl;
		}
    }
    
    /*  Creates a map to keep track of how much of each vertex weight is 
        currently being used by incident edges. Each pair's first element is 
        the currently used weight, and the second element is the total vertex 
        weight */
    map<vertex_id, pair<unsigned int, unsigned int> > used_weights;
    
    // Initialize the used_weights map
    for (map<vertex_id, unsigned int>::const_iterator weight_it = weights.begin();
        weight_it != weights.end(); weight_it++)
    {
        vertex_id current_vertex = weight_it->first;
        unsigned int vertex_weight = weight_it->second;
        used_weights[current_vertex] = make_pair(0, vertex_weight);
    }

    for (vector<edge>::const_reverse_iterator edge_it = edge_vec.rbegin();
        edge_it != edge_vec.rend(); edge_it++)
    {
        vertex_id v1 = edge_it->first;
        vertex_id v2 = edge_it->second;

        pair<unsigned int, unsigned int>& used_weight_v1 = used_weights[v1];
        pair<unsigned int, unsigned int>& used_weight_v2 = used_weights[v2];

        // Compute the remaining price that edges incident in v1 and v2 can have
        unsigned int available_price_v1 = used_weight_v1.second - used_weight_v1.first;
        unsigned int available_price_v2 = used_weight_v2.second - used_weight_v2.first;

		if (DEBUG_MESSAGES)
		{
			cout << "Edge: (" << v1 << ", " << v2 << ")" << endl;
			cout << "    Price sum (" << v1 << "): " << used_weight_v1.first << "/" << used_weight_v1.second << endl;
			cout << "    Price sum (" << v2 << "): " << used_weight_v2.first << "/" << used_weight_v2.second << endl;
			cout << "    Available price (" << v1 << "): " << available_price_v1 << endl;
			cout << "    Available price (" << v2 << "): " << available_price_v2 << endl;
		}

        // Ignore this edge if either v1 or v2 is tight
        if ( available_price_v1 == 0 || available_price_v2 == 0 )
        {
			if (DEBUG_MESSAGES)
				cout << "    **** Tight vertex ****" << endl;
            continue;
        }

        /*  Choose the price of this edge to be the minimum available price, so
            that the fairness condition is preserved */
        unsigned int price = min( available_price_v1, available_price_v2 );

        // Increase the tightness of v1 and v2 by the price of the edge
        used_weight_v1.first += price;
        used_weight_v2.first += price;

		if (DEBUG_MESSAGES)
		{
			cout << "    Edge price: " << price << endl;
			cout << "    New price sum (" << v1 << "): " << used_weight_v1.first << "/" << used_weight_v1.second << endl;
			cout << "    New price sum (" << v2 << "): " << used_weight_v2.first << "/" << used_weight_v2.second << endl;
		}
    }

    cover.clear();
    
	// Counter for the total weight sum of the vertex cover
	unsigned int cover_weight_sum = 0;
	
	// Add all tight vertices to the vertex cover
    for (auto used_weights_it = used_weights.begin(); 
        used_weights_it != used_weights.end(); used_weights_it++)
    {
        vertex_id current_vertex = used_weights_it->first;
        unsigned int used_weight = used_weights_it->second.first;
        unsigned int vertex_weight = used_weights_it->second.second;

        // If the vertex is tight, add it to the vertex cover
		if (used_weight == vertex_weight)
		{
            cover.insert( current_vertex );
			cover_weight_sum += vertex_weight;
		}
    }

	return cover_weight_sum;
}

// ============================================================================



// ============================================================================
// =========================== ANALYSIS FUNCTIONS =============================
// ============================================================================

struct Statistics {
	double sum_greedy, sum_pricing, sum_improved_pricing;
	double worst_greedy_dist, worst_pricing_dist, worst_improved_pricing_dist;

	Statistics()
	{
		sum_greedy = sum_pricing = sum_improved_pricing = 0;
		worst_greedy_dist = worst_pricing_dist = worst_improved_pricing_dist = 0;
	}
};

template <class T>
string set_to_string(const T& set)
{
	stringstream ss;
	
	ss << "{";
	for (T::const_iterator it = set.begin(); it != set.end(); it++)
	{
		if (it != set.begin()) { ss << ", "; }
		ss << *it;
	}
	ss << "}";
	
	return ss.str();
}

void analyse_input(ostream* output, const Simple_Graph& graph)
{
	set<vertex_id> optimal_cover, greedy_cover, pricing_cover, improved_pricing_cover;
	unsigned int optimal_weight_sum, greedy_weight_sum, pricing_weight_sum, improved_pricing_weight_sum;

	optimal_weight_sum = graph.optimal_vertex_cover(optimal_cover);
	greedy_weight_sum = graph.greedy_vertex_cover(greedy_cover);
	pricing_weight_sum = graph.pricing_vertex_cover(pricing_cover);
	improved_pricing_weight_sum = graph.improved_pricing_vertex_cover(improved_pricing_cover);

	stringstream ss;

	// Print results for the optimal vertex cover
	for (set<vertex_id>::const_iterator it = optimal_cover.begin(); it != optimal_cover.end(); it++)
		ss << *it << " ";
	ss << optimal_weight_sum << endl;

	// Print results for the greedy vertex cover
	for (set<vertex_id>::const_iterator it = greedy_cover.begin(); it != greedy_cover.end(); it++)
		ss << *it << " ";
	ss << greedy_weight_sum << endl;

	// Print results for the pricing vertex cover
	for (set<vertex_id>::const_iterator it = pricing_cover.begin(); it != pricing_cover.end(); it++)
		ss << *it << " ";
	ss << pricing_weight_sum << endl;

	// Print results for the improved_pricing vertex cover
	for (set<vertex_id>::const_iterator it = improved_pricing_cover.begin(); it != improved_pricing_cover.end(); it++)
		ss << *it << " ";
	ss << improved_pricing_weight_sum << endl;

	ss << endl;

	ostream_lock.lock();
	*output << ss.str();
	ostream_lock.unlock();
}

void analyse_simulation(ostream* os, Statistics* statistics, const Simple_Graph graph)
{
	set<vertex_id> optimal_cover, greedy_cover, pricing_cover, improved_pricing_cover;
	double optimal_weight_sum, greedy_weight_sum, pricing_weight_sum, improved_pricing_weight_sum;

	optimal_weight_sum = graph.optimal_vertex_cover(optimal_cover);
	greedy_weight_sum = graph.greedy_vertex_cover(greedy_cover);
	pricing_weight_sum = graph.pricing_vertex_cover(pricing_cover);
	improved_pricing_weight_sum = graph.improved_pricing_vertex_cover(improved_pricing_cover);

	// Compute distances to optimal
	double dist_greedy_optimal = greedy_weight_sum / optimal_weight_sum;
	double dist_pricing_optimal = pricing_weight_sum / optimal_weight_sum;
	double dist_improved_pricing_optimal = improved_pricing_weight_sum / optimal_weight_sum;

	// Update counters for average
	statistics->sum_greedy += dist_greedy_optimal;
	statistics->sum_pricing += dist_pricing_optimal;
	statistics->sum_improved_pricing += dist_improved_pricing_optimal;

	// Update worst distances
	statistics->worst_greedy_dist = max( statistics->worst_greedy_dist, dist_greedy_optimal );
	statistics->worst_pricing_dist = max( statistics->worst_pricing_dist, dist_pricing_optimal );
	statistics->worst_improved_pricing_dist = max( statistics->worst_improved_pricing_dist, 
		dist_improved_pricing_optimal );

	ostream_lock.lock();

	// Print the graph to the output stream
	*os << "******************************" << endl << endl;
	*os << "=== Graph ===" << endl << graph.to_string_output() << endl;
	// Print the solutions to the output stream
	*os << "=== Solutions ===" << endl;
	*os << "Optimal: " << set_to_string(optimal_cover) << endl;
	*os << "Greedy: " << set_to_string(greedy_cover) << endl;
	*os << "Pricing: " << set_to_string(pricing_cover) << endl;
	*os << "Improved Pricing: " << set_to_string(improved_pricing_cover) << endl << endl;
	// Print the distances from the heuristic solutions to OPT
	*os << "=== Distances to OPT ===" << endl;
	*os << "Greedy: " << dist_greedy_optimal << endl;
	*os << "Pricing: " << dist_pricing_optimal << endl;
	*os << "Improved Pricing: " << dist_improved_pricing_optimal << endl << endl;

	ostream_lock.unlock();
}

void run_simulation(ostream& os, unsigned int num_vertices, 
	unsigned int num_edges, unsigned int max_weight)
{
	Statistics statistics;
	stringstream thread_ostream;

	vector<thread> running_threads;
	for (int i = 0; i < NUM_SIMULATIONS; i++)
	{
		Simple_Graph graph;
		graph.generate_as_random(num_vertices, num_edges, max_weight);

		running_threads.push_back( thread(analyse_simulation, &thread_ostream, 
			&statistics, graph) );
		
		// If NUM_THREADS threads are running, wait for them and then clear the vector
		if (running_threads.size() == NUM_THREADS)
		{
			for (vector<thread>::size_type i = 0; i < running_threads.size(); i++)
				running_threads[i].join();

			running_threads.clear();
		}
	}

	// Wait for remaining threads
	for (vector<thread>::size_type i = 0; i < running_threads.size(); i++)
		running_threads[i].join();

	// Print simulation statistics to the output stream
	os << "========== GLOBAL STATISTICS ==========" << endl;
	os << "Average greedy distance to OPT: " << statistics.sum_greedy / NUM_SIMULATIONS << endl;
	os << "Worst greedy distance to OPT: " << statistics.worst_greedy_dist << endl << endl;
	os << "Average pricing distance to OPT: " << statistics.sum_pricing / NUM_SIMULATIONS << endl;
	os << "Worst pricing distance to OPT: " << statistics.worst_pricing_dist << endl << endl;
	os << "Average improved pricing distance to OPT: " << statistics.sum_improved_pricing / NUM_SIMULATIONS << endl;
	os << "Worst improved pricing distance to OPT: " << statistics.worst_improved_pricing_dist << endl;
	os << "=======================================" << endl << endl;

	// Print thread outputs to the output stream
	os << thread_ostream.str() << endl;
}

void parse_input(istream& is, Simple_Graph& graph)
{
	graph.clear();

	unsigned int vertex_num;
	
	if ( !(is >> vertex_num) ) throw runtime_error("Could not read from stream");
	
	for (unsigned int vertex = 1; vertex <= vertex_num; vertex++)
		graph.add_vertex(vertex, 1);

	for (unsigned int vertex = 1; vertex <= vertex_num; vertex++)
	{
		unsigned int weight;
		if ( !(is >> weight) ) throw runtime_error("Could not read from stream");

		graph.set_weight(vertex, weight);

		unsigned int adjacent;
		while ( is >> adjacent && adjacent != 0 )
			graph.add_edge(vertex, adjacent);
	}
}

void process_input(istream& is, ostream& os)
{
	unsigned int num_inputs;
	if ( !(is >> num_inputs) ) throw runtime_error("Could not read from file");

	vector<thread> running_threads;
	
	for (unsigned int instance = 0; instance < num_inputs; instance++)
	{
		Simple_Graph graph;
		parse_input(is, graph);

		running_threads.push_back( thread(analyse_input, &os, graph) );

		// If NUM_THREADS threads are running, wait for them and then clear the vector
		if (running_threads.size() == NUM_THREADS)
		{
			for (vector<thread>::size_type i = 0; i < running_threads.size(); i++)
				running_threads[i].join();

			running_threads.clear();
		}
	}

	// Wait for remaining threads
	for (vector<thread>::size_type i = 0; i < running_threads.size(); i++)
		running_threads[i].join();
}

// ============================================================================



// ============================================================================
// ============================== MAIN FUNCTION ===============================
// ============================================================================


int main()
{
	try
	{
		string mode;
		cout << "Type the program mode: \"simulation\" or \"input\"..." << endl;
		while ( getline(cin, mode) && mode != "simulation" && mode != "input" )
		{
			cout << "Invalid mode" << endl;
		}

		if ( !cin )
			throw runtime_error("Cannot read from standard input stream");

		stringstream analysis_output;

		if (mode == "input")
		{
			stringstream input_file_content;

			// ===== Read input file =====
			ifstream input_file("WVC_jgsma.txt");

			string line;
			while (getline(input_file, line))
				input_file_content << line << endl;

			input_file.close();
			// ===========================

			// Process input
			process_input(input_file_content, analysis_output);
		}
		else
		{
			// ========== Read the simulation parameters ==========
			int num_vertices, num_edges, max_weight;
			cout << "Type the number of vertices in each graph..." << endl;
			while (!(cin >> num_vertices) || num_vertices <= 0)
			{
				cout << "Invalid number of vertices" << endl;
				cin.clear();
			}

			cout << "Type the number of edges in each graph..." << endl;
			while (!(cin >> num_edges) || num_edges <= 0)
			{
				cout << "Invalid number of edges" << endl;
				cin.clear();
			}

			cout << "Type the max weight for each vertex..." << endl;
			while (!(cin >> max_weight) || max_weight <= 0)
			{
				cout << "Invalid weight value" << endl;
				cin.clear();
			}

			cout << "Type the number of instances in the simulation..." << endl;
			while (!(cin >> NUM_SIMULATIONS) || NUM_SIMULATIONS <= 0)
			{
				cout << "Invalid number of instances" << endl;
				cin.clear();
			}
			// ====================================================

			run_simulation(analysis_output, num_vertices, num_edges, max_weight);
		}

		// ===== Write results to output file =====
		ofstream output_file("WVC_jgsma_Output.txt");
		output_file << analysis_output.str() << endl;
		output_file.close();
		// ========================================

	}
	catch (std::exception& e)
	{
		cout << "ERROR: " << e.what() << endl;
		return 1;
	}

	return 0;
}
// ============================================================================
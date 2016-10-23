#include <algorithm>
#include <fstream>
#include <iostream>
#include <map>
#include <random>
#include <set>
#include <stdexcept>
#include <sstream>
#include <utility>
#include <vector>

using std::cout;
using std::endl;
using std::fstream;
using std::logic_error;
using std::make_pair;
using std::map;
using std::max;
using std::min;
using std::mt19937;
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
using std::uniform_int_distribution;
using std::vector;

typedef unsigned int vertex_id;

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

    void greedy_vertex_cover(set<vertex_id> &cover) const;

    void pricing_vertex_cover(set<vertex_id>& cover) const;

    void improved_pricing_vertex_cover(set<vertex_id>& cover) const;

    string to_string() const;

	string to_string_output() const;

    void generate_as_random(unsigned int num_vertices, unsigned int num_edges,
        unsigned int max_weight);


private:
    map<vertex_id, set<vertex_id> > graph;
    map<vertex_id, unsigned int> weights;

    bool equal_weights(const pair<vertex_id, vertex_id>& edge) const;

    vertex_id greatest_ratio_vertex(
        const map<vertex_id, set<vertex_id> > &graph) const;

    void pricing_vertex_cover( set<vertex_id>& cover, bool improved ) const;

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
        if ( i != graph.begin() ) { ss << endl; };

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
        
    }

    return ss.str();
}


string Simple_Graph::to_string_output() const
{
	stringstream ss;
	for (map<vertex_id, set<vertex_id> >::const_iterator i = graph.begin();
		i != graph.end(); i++)
	{
		if (i != graph.begin()) { ss << endl; };

		vertex_id current = i->first;
		const set<vertex_id>& adjacent = i->second;

		for (set<vertex_id>::const_iterator j = adjacent.begin();
			j != adjacent.end(); j++)
		{
			if (j != adjacent.begin())
				ss << " ";

			ss << *j;
		}

	}

	return ss.str();
}


vertex_id Simple_Graph::greatest_ratio_vertex(
    const map<vertex_id, set<vertex_id> > &graph) const
{
    double best_ratio = -1;
    vertex_id best_vertex;
    typedef map<vertex_id, set<vertex_id> >::const_iterator vertex_iterator;
    
    #ifdef _DEBUG
    cout << "Vertices:" << endl;
    #endif

    for(vertex_iterator v_it = graph.begin(); v_it != graph.end(); v_it++)
    {
        vertex_id current_vertex = v_it->first;
        set<vertex_id>::size_type degree = v_it->second.size();

		double vertex_weight = weights.find(current_vertex)->second;
		double ratio = degree / vertex_weight;

        #ifdef _DEBUG
        cout << "    Vertex: " << current_vertex << endl;
        cout << "    Degree: " << degree << endl;
        cout << "    Weight: " << vertex_weight << endl;
        cout << "    Ratio: " << ratio << endl;
        cout << endl;
        #endif

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

unsigned int rand_uint(unsigned int range_limit)
{
	uniform_int_distribution<unsigned int> dist(0, range_limit - 1);
	return dist(mt);
}


void Simple_Graph::greedy_vertex_cover(set<vertex_id>& cover) const
{
    cover.clear();

    map<vertex_id, set<vertex_id> > graph_copy = graph;

    while( !graph_copy.empty() )
    {
        // Find the vertex with the greatest degree/weight ratio
        vertex_id grv = greatest_ratio_vertex(graph_copy);

        #ifdef _DEBUG
        cout << "Greatest ratio vertex: " << grv << endl;
        cout << "********************" << endl << endl;
        #endif

        cover.insert( grv );
        
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
}


void Simple_Graph::pricing_vertex_cover(set<vertex_id>& cover) const
{
    pricing_vertex_cover(cover, false); 
}


void Simple_Graph::improved_pricing_vertex_cover(set<vertex_id>& cover) const
{
    pricing_vertex_cover(cover, true); 
}


void Simple_Graph::pricing_vertex_cover( set<vertex_id>& cover, bool improved) const
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

#ifdef _DEBUG
    cout << "Edge vector:" << endl;
    for (auto it = edge_vec.begin(); it != edge_vec.end(); it++)
        cout << "  (" << it->first << ", " << it->second << ")" << endl;
#endif

    random_shuffle(edge_vec.begin(), edge_vec.end(), rand_uint);

#ifdef _DEBUG
    cout << "Shuffled edge vector:" << endl;
    for (auto it = edge_vec.begin(); it != edge_vec.end(); it++)
        cout << "  (" << it->first << ", " << it->second << ")" << endl;
#endif

    /*  Improved pricing heuristic. This heuristic partitions the vector in such
        a way that all edges whose vertices have different weights will be chosen
        before edges whose vertices have equal weights */
    if (improved)
    {
        stable_partition(edge_vec.begin(), edge_vec.end(), 
            [this] (const pair<vertex_id, vertex_id>& edge) { return this->equal_weights(edge); } );
        
#ifdef _DEBUG
        cout << "Edge vector after heuristic:" << endl;
        for (auto it = edge_vec.begin(); it != edge_vec.end(); it++)
            cout << "  (" << it->first << ", " << it->second << ")" << endl;
#endif
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

#ifdef _DEBUG
        cout << "Edge: (" << v1 << ", " << v2 << ")" << endl;
        cout << "    Price sum (" << v1 << "): " << used_weight_v1.first << "/" << used_weight_v1.second << endl;
        cout << "    Price sum (" << v2 << "): " << used_weight_v2.first << "/" << used_weight_v2.second << endl;
        cout << "    Available price (" << v1 << "): " << available_price_v1 << endl;
        cout << "    Available price (" << v2 << "): " << available_price_v2 << endl;
#endif

        // Ignore this edge if either v1 or v2 is tight
        if ( available_price_v1 == 0 || available_price_v2 == 0 )
        {
#ifdef _DEBUG
            cout << "    **** Tight vertex ****" << endl;
#endif
            continue;
        }

        /*  Choose the price of this edge to be the minimum available price, so
            that the fairness condition is preserved */
        unsigned int price = min( available_price_v1, available_price_v2 );

        // Increase the tightness of v1 and v2 by the price of the edge
        used_weight_v1.first += price;
        used_weight_v2.first += price;

#ifdef _DEBUG
        cout << "    Edge price: " << price << endl;
        cout << "    New price sum (" << v1 << "): " << used_weight_v1.first << "/" << used_weight_v1.second << endl;
        cout << "    New price sum (" << v2 << "): " << used_weight_v2.first << "/" << used_weight_v2.second << endl;
#endif
    }

    cover.clear();
    
    for (auto used_weights_it = used_weights.begin(); 
        used_weights_it != used_weights.end(); used_weights_it++)
    {
        vertex_id current_vertex = used_weights_it->first;
        unsigned int used_weight = used_weights_it->second.first;
        unsigned int vertex_weight = used_weights_it->second.second;

        // If the vertex is tight, add it to the vertex cover
        if ( used_weight == vertex_weight )
            cover.insert( current_vertex );
    }
}

// ============================================================================



// ============================================================================
// =========================== ANALYSIS FUNCTIONS =============================
// ============================================================================

// ============================================================================



// ============================================================================
// ============================== MAIN FUNCTION ===============================
// ============================================================================


int main()
{
	Simple_Graph g;

    g.generate_as_random(5, 8, 10);

    cout << g.to_string() << endl;
        
    set<vertex_id> cover;
    g.improved_pricing_vertex_cover(cover);
    
    cout << "Vertex cover:" << endl;
    for (auto it = cover.begin(); it != cover.end(); it++)
    {
        cout << *it << " ";
    }
    cout << endl;

    return 0;
}
// ============================================================================
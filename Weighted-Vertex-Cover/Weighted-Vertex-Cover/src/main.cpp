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

private:
    map<vertex_id, set<vertex_id> > graph;
    map<vertex_id, unsigned int> weights;

    bool is_valid_cover(const set<vertex_id> &combination) const;
    
    vertex_id greatest_ratio_vertex(
        const map<vertex_id, set<vertex_id> > &graph) const;
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
        const set<vertex_id>& adjacent = i->second;

        ss << current << ": ";

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
    double best_ratio = 0;
    vertex_id best_vertex;
    typedef map<vertex_id, set<vertex_id> >::const_iterator vertex_iterator;
        
    for(vertex_iterator v_it = graph.begin(); v_it != graph.end(); v_it++)
    {
        vertex_id current_vertex = v_it->first;
        set<vertex_id>::size_type degree = v_it->second.size();

		unsigned int vertex_weight = weights.find(current_vertex)->second;
		double ratio = degree / vertex_weight;

        if ( ratio > best_ratio )
        {
            best_ratio = ratio;
            best_vertex = current_vertex;
        }
    }

    return best_vertex;
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
        vertex_id gdv = greatest_ratio_vertex(graph_copy);
        
        cover.insert( gdv );
        
        /*  Remove the greatest ratio vertex from the adjacency list of each
            of its adjacent vertices */
        set<vertex_id> &adj_vertices = graph_copy[gdv];
        for (set<vertex_id>::const_iterator adj_vertex_it = adj_vertices.begin();
            adj_vertex_it != adj_vertices.end(); adj_vertex_it++)
        {
            vertex_id adj_vertex = *adj_vertex_it;
            set<vertex_id> &adj_vertex_edges = graph_copy[adj_vertex];
            
            adj_vertex_edges.erase( gdv );
            
            // If the vertex has become isolated, remove it from graph_copy
            if ( adj_vertex_edges.empty() )
                graph_copy.erase( *adj_vertex_it );
        }

        // Remove the greatest degree vertex from the graph
        graph_copy.erase( gdv );
    }
}


void Simple_Graph::pricing_vertex_cover(set<vertex_id>& cover) const
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

    random_shuffle(edge_vec.begin(), edge_vec.end(), rand_uint);
    
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
 
        // Ignore this edge if either v1 or v2 is tight
        if ( available_price_v1 == 0 || available_price_v2 == 0 )
            continue;

        /*  Choose the price of this edge to be the minimum available price, so
            that the fairness condition is preserved */
        unsigned int price = min( available_price_v1, available_price_v2 );

        // Increase the tightness of v1 and v2 by the price of the edge
        used_weight_v1.first += price;
        used_weight_v2.first += price;
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


void Simple_Graph::improved_pricing_vertex_cover(set<vertex_id>& cover) const
{
    // TODO...
}


bool Simple_Graph::is_valid_cover(const set<vertex_id> &combination) const
{
    for (map<vertex_id, set<vertex_id> >::const_iterator it = graph.begin();
        it != graph.end(); it++)
    {
        vertex_id current_vertex = it->first;
        const set<vertex_id> &edges = it->second;

        // Continue if the current vertex is in combination
        if ( combination.find(current_vertex) != combination.end() )
            continue;

        for (set<vertex_id>::const_iterator edge_it = edges.begin();
            edge_it != edges.end(); edge_it++)
        {
            vertex_id other_vertex = *edge_it;

            // Check if the second vertex is in combination
            if ( combination.find(other_vertex) == combination.end() )
                return false;   // Found an edge not covered           
        }
    }

    return true;
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
	return 0;
}
// ============================================================================
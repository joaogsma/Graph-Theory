#include <cstdlib>
#include <iostream>
#include <map>
#include <set>
#include <stdexcept>
#include <sstream>
#include <time.h>
#include <vector>

using std::cout;
using std::endl;
using std::logic_error;
using std::map;
using std::string;
using std::stringstream;
using std::set;
using std::vector;


// ============================================================================
// ============================= CLASS DEFINITION =============================
// ============================================================================

typedef unsigned int vertex_id;

class Simple_Graph {
public:
    void clear() { graph.clear(); }

    void clear_edges();
    
    void add_vertex(vertex_id id);

    void add_edge(vertex_id a, vertex_id b);

    void greedy_vc(set<vertex_id> &cover);

    void optimal_vc(set<vertex_id> &cover);

    string to_string();

private:
    map<vertex_id, set<vertex_id> > graph;

    bool is_valid_cover(const set<vertex_id> &combination);
    
    vertex_id greatest_degree_vertex(
        const map<vertex_id, set<vertex_id> > &graph);
};

// ============================================================================



// ============================================================================
// =========================== FUNCTION DEFINITIONS ===========================
// ============================================================================

string Simple_Graph::to_string()
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


void Simple_Graph::clear_edges()
{
    for(map<vertex_id, set<vertex_id> >::iterator it = graph.begin();
        it != graph.end(); it++)
    {
        it->second.clear();
    } 
}


void Simple_Graph::add_vertex(vertex_id id)
{
    // Check if the vertex has already been added
    if ( graph.find(id) != graph.end() )
    {
        stringstream ss;
        ss << "Vertex " << id << " already in the graph" << endl;
        throw logic_error( ss.str() );
    }

    // Create the vertex in the graph, with no edges
    graph[id];
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


vertex_id Simple_Graph::greatest_degree_vertex(
    const map<vertex_id, set<vertex_id> > &graph)
{
    set<vertex_id>::size_type max_degree = 0;
    vertex_id greatest_deg_vertex;
    typedef map<vertex_id, set<vertex_id> >::const_iterator vertex_iterator;
        
    for(vertex_iterator v_it = graph.begin(); 
        v_it != graph.end(); v_it++)
    {
        vertex_id current_vertex = v_it->first;
        set<vertex_id>::size_type degree = v_it->second.size();
        
        if ( degree > max_degree )
        {
            max_degree = degree;
            greatest_deg_vertex = current_vertex;
        }
    }

    return greatest_deg_vertex;
}


void Simple_Graph::greedy_vc(set<vertex_id> &cover)
{
    cover.clear();

    map<vertex_id, set<vertex_id> > graph_copy = graph;

    while( !graph_copy.empty() )
    {
        // Find the vertex with the greatest degree
        vertex_id gdv = greatest_degree_vertex(graph_copy);
        
        cover.insert( gdv );
        
        /*  Remove the greatest degree vertex from the adjacency list of each
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


void Simple_Graph::optimal_vc(set<vertex_id> &cover)
{
    cover.clear();

    // This set will hold the vertex sets checked as possible covers
    set<set<vertex_id> > valid_covers;
    set<set<vertex_id> > *subsets = new set<set<vertex_id> >;
    set<set<vertex_id> > *next_subsets = new set<set<vertex_id> >;

    // Insert all sets of size 1 in subsets
    for (map<vertex_id, set<vertex_id> >::const_iterator it = graph.begin();
        it != graph.end(); it++)
    {
        set<vertex_id> unit_set;
        unit_set.insert( it->first );

        // If the unit set is a valid cover, return it
        if ( is_valid_cover(unit_set) )
        {
            cover = unit_set;
            delete subsets;
            delete next_subsets;
            return;
        }

        subsets->insert( unit_set );
    }

    while(true)
    {
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
                if ( current_subset.find( current_vertex ) != current_subset.end() )
                    continue;

                set<vertex_id> current_subset_copy = current_subset;

                /*  Add the current subset U {current_vertex} to the next
                    iteration set */
                current_subset_copy.insert( current_vertex );

                // Check if the subset is a valid vertex cover
                if ( is_valid_cover(current_subset_copy) )
                {
                    cover = current_subset_copy;
                    delete subsets;
                    delete next_subsets;
                    return;
                }

                next_subsets->insert( current_subset_copy );
            }
        }
        
        // Update the subsets pointer
        delete subsets;
        subsets = next_subsets;
        next_subsets = new set<set<vertex_id> >; 
    }
}


bool Simple_Graph::is_valid_cover(const set<vertex_id> &combination)
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
// ======================== GRAPH GENERATION FUNCTIONS ========================
// ============================================================================

vertex_id pick_random(const set<vertex_id>& vertex_set)
{
    const set<vertex_id>::size_type rand_idx = rand() % vertex_set.size();

    // Advance it to the element at position rand_idx
    set<vertex_id>::const_iterator it = vertex_set.begin();
    for ( set<vertex_id>::size_type i = 0; i < rand_idx; i++ )
        it++;

    return *it;
}


vector<vertex_id> loop_erasing_random_walk(
    const set<vertex_id>& tree_vertices, const set<vertex_id>& vertex_set)
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


void wilson_spanning_tree(Simple_Graph& graph, unsigned int num_vertices)
{
    graph.clear_edges();

    // Create the vertex set
    set<vertex_id> vertex_set;
    for (unsigned int i = 1; i <= num_vertices; i++)
        vertex_set.insert(i);

    // Create the set of vertices used in the spanning tree
    set<vertex_id> tree_vertices;
    // Add vertex 1 as the root
    tree_vertices.insert(1);

    vector<vertex_id> walk = loop_erasing_random_walk(tree_vertices, vertex_set);
    while( !walk.empty() )
    {
        tree_vertices.insert( walk.begin(), walk.end() );

        for (vector<vertex_id>::size_type i = 0; i < walk.size() - 1; i++)
            graph.add_edge( walk[i], walk[i+1] );

        walk = loop_erasing_random_walk(tree_vertices, vertex_set);
    }
}


void add_random_edges(Simple_Graph& graph, unsigned int num_vertices, 
    unsigned int target_num_edges)
{

}


void generate_graph(Simple_Graph& graph, unsigned int num_vertices, 
    unsigned int num_edges)
{
    if (num_edges < num_vertices - 1)
    {
        stringstream ss;
        ss << "A graph with |V| = " << num_vertices << " and |E| = " << 
            num_edges << "cannot be connected" << endl;
        throw logic_error( ss.str() );
    }

    // Add the specified number of vertices
    for (unsigned int i = 1; i <= num_vertices; i++)
        graph.add_vertex(i);

    // Generate a random spanning tree to ensure the graph is connected
    wilson_spanning_tree(graph, num_vertices);

    // Add random edges to achieve the specified density
    add_random_edges(graph, num_vertices, num_edges);
}

// ============================================================================



// ============================================================================
// ============================== MAIN FUNCTION ===============================
// ============================================================================

int main()
{

    int n = 10;

    srand( time(NULL) );

    Simple_Graph g;
    for (vertex_id i = 1; i <= n; i++)
        g.add_vertex(i);

    wilson_spanning_tree(g, n);

    cout << g.to_string();

    return 0;
}

// ============================================================================
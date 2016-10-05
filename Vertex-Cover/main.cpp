#include <iostream>
#include <map>
#include <set>
#include <stdexcept>
#include <sstream>

using std::cout;
using std::endl;
using std::logic_error;
using std::map;
using std::string;
using std::stringstream;
using std::set;


// ============================================================================
// ============================= CLASS DEFINITION =============================
// ============================================================================

typedef unsigned int vertex_id;

class Simple_Graph {
private:
    map<vertex_id, set<vertex_id> > graph;

    bool is_valid_cover(const set<vertex_id> &combination);
    
    vertex_id greatest_degree_vertex(
        const map<vertex_id, set<vertex_id> > &graph);

public:
    void add_vertex(vertex_id id);

    void add_edge(vertex_id a, vertex_id b);

    void greedy_vc(set<vertex_id> &cover);

    void optimal_vc(set<vertex_id> &cover);

};

// ============================================================================



// ============================================================================
// =========================== FUNCTION DEFINITIONS ===========================
// ============================================================================

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

// TODO...

// ============================================================================



// ============================================================================
// ============================== MAIN FUNCTION ===============================
// ============================================================================

int main()
{
    // Test sample
    Simple_Graph g;

    int n = 30;

    for (int i = 1; i <= n; i++)
        g.add_vertex(i);

    for (int i = 1; i <= n; i++)
    {
        for (int j = 1; j <= n; j++)
            g.add_edge(i, j);
    }

    set<vertex_id> cover;

    g.optimal_vc( cover );

    cout << "Vertex Cover: {";
    
    set<vertex_id>::const_iterator it = cover.begin();
    cout << *it;
    it++;

    while( it != cover.end() )
    {
        cout << ", " << *it;
        it++;
    }
    cout << "}" << endl;

    return 0;
}

// ============================================================================
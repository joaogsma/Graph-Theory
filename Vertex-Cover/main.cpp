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

/*

===== TODO List: =====
- Greedy vertex cover function (heuristic)
- Main function
- Graph generation functions
- Analysis functions
======================

*/


// ============================================================================
// ============================= CLASS DEFINITION =============================
// ============================================================================

typedef unsigned int vertex_id;

class Simple_Graph {
private:
    map<vertex_id, set<vertex_id> > graph;

public:
    void add_vertex(vertex_id id);

    void add_edge(vertex_id a, vertex_id b);

    void greedy_vc(set<vertex_id> &cover);

    void optimal_vc(set<vertex_id> &cover);

    bool is_valid_cover(const set<vertex_id> &combination);
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


void Simple_Graph::greedy_vc(set<vertex_id> &cover)
{
    cover.clear();
}


void Simple_Graph::optimal_vc(set<vertex_id> &cover)
{
    cover.clear();

    // This set will hold the vertex sets checked as possible covers
    set<set<vertex_id> > valid_covers;
    set<set<vertex_id> > *combinations = new set<set<vertex_id> >;
    set<set<vertex_id> > *next_combinations = new set<set<vertex_id> >;

    // Insert all sets of size 1 in combinations
    for (map<vertex_id, set<vertex_id> >::const_iterator it = graph.begin();
        it != graph.end(); it++)
    {
        set<vertex_id> temp_set;
        temp_set.insert( it->first );
        combinations->insert( temp_set );
    }

    while(true)
    {
        cout << "Combinations:" << endl;
        for (set<set<vertex_id> >::const_iterator it = combinations->begin();
            it != combinations->end(); it++)
        {
            cout << "\t";
            for (set<vertex_id>::const_iterator id_it = it->begin(); 
                id_it != it->end(); id_it++)
            {
                cout << " " << *id_it;
            }
            cout << endl;
        }

        for (set<set<vertex_id> >::const_iterator comb_it = combinations->begin();
            comb_it != combinations->end(); comb_it++)
        {
            const set<vertex_id> &current_comb = *comb_it;

            // For each vertex and each subset, add a subset that contains the vertex 
            // and one that does not
            for (map<vertex_id, set<vertex_id> >::const_iterator vertex_it = graph.begin();
                vertex_it != graph.end(); vertex_it++)
            {
                vertex_id current_vertex = vertex_it->first;
           
                // Continue if the current vertex is already in it
                if ( current_comb.find( current_vertex ) != current_comb.end() )
                    continue;

                set<vertex_id> current_comb_copy = current_comb;

                // Add the current combination U {current_vertex} to the next
                // iteration set
                current_comb_copy.insert( current_vertex );

                // Check if the combination is a valid vertex cover
                if ( is_valid_cover(current_comb_copy) )
                    cover = current_comb_copy;

                next_combinations->insert( current_comb_copy );
            }

            // Stop if a vertex cover was found
            if ( !cover.empty() )
                return;
        }
        
        // update the combinations pointer
        delete combinations;
        combinations = next_combinations;
        next_combinations = new set<set<vertex_id> >; 
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

    for (int i = 1; i <= 6; i++)
        g.add_vertex(i);

    g.add_edge(1, 2);
    g.add_edge(1, 3);
    g.add_edge(2, 3);
    g.add_edge(2, 4);
    g.add_edge(2, 5);
    g.add_edge(2, 6);
    
    set<vertex_id> cover;

    g.optimal_vc( cover );

    cout << "Vertex Cover:" << endl;
    for (set<vertex_id>::const_iterator it = cover.begin(); it != cover.end(); it++)
    {
        cout << " " << *it;
    }
    cout << endl;

    return 0;
}

// ============================================================================
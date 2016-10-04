#include <iostream>
#include <map>
#include <set>
#include <stdexcept>
#include <sstream>

using std::endl;
using std::logic_error;
using std::map;
using std::string;
using std::stringstream;
using std::set;
using std::cout;

/*

===== TODO List: =====
- Optimal vertex cover function (brute-force)
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
    Simple_Graph();

    void add_vertex(vertex_id id);

    void add_edge(vertex_id a, vertex_id b);

    void optimal_vc(set<vertex_id> &cover);

    void greedy_vc(set<vertex_id> &cover);
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


void Simple_Graph::optimal_vc(set<vertex_id> &cover)
{
    cover.clear();

    // This set will hold the vertex sets checked as possible covers
    set<set<vertex_id> > valid_covers;
    set<set<vertex_id> > combinations;

    // Insert all sets of size 1 in combinations
    for (map<vertex_id, set<vertex_id> >::const_iterator it = graph.begin();
        it != graph.end(); it++)
    {
        set<vertex_id> temp_set;
        temp_set.insert( it->first );
        combinations.insert( temp_set );
    }

    // For each vertex and each subset, add a subset that contains the vertex 
    // and one that does not
    for (map<vertex_id, set<vertex_id> >::const_iterator it = graph.begin();
        it != graph.end(); it++)
    {
        for (set<set<vertex_id> >::const_iterator)
    }
}





void Simple_Graph::greedy_vc(set<vertex_id> &cover)
{
    cover.clear();


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
    cout << "Hello, World!" << endl;

    return 0;
}

// ============================================================================
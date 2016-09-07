#include <exception>
#include <fstream>
#include <iostream>
#include <map>
#include <stdexcept>
#include <string>
#include <sstream>
#include <utility>
#include <vector>

using std::cout;            using std::logic_error;
using std::domain_error;    using std::map;
using std::endl;            using std::pair;
using std::exception;       using std::runtime_error;
using std::fstream;         using std::stringstream;
using std::getline;         using std::string;
using std::istream;         using std::vector;
using std::logic_error;

struct Edge {
    int flow, capacity;

    Edge() : flow(0), capacity(0) {}
    Edge(int flow, int capacity) : flow(flow), capacity(capacity) {}
};

class Network {
private:
    map< unsigned int, map<unsigned int, Edge> > net;

    
    void residual_network(Network& res_net) {
        // TODO
    }


    void check_edge(unsigned int orig, unsigned int dest) {
        bool exists = true;
        
        // Check weather vertex orig exists
        if (net.find(orig) == net.end())
            exists = false;
        
        // Check weather vertex dest exists
        if (net.find(dest) == net.end())
            exists = false;
        
        // Check weather there is an edge between them
        if (net[orig].find(dest) == net[orig].end())
            exists = false;

        // If edge does not exist, throw an exception
        if (!exists) {
            stringstream ss;
            ss << "There is no edge (" << orig << "," << dest << ")";
            throw logic_error( ss.str() );
        }

    }

public:
    int source, sink;


    // ======================= MISCELLANEOUS FUNCTIONS ========================

    void clean_flow() {
        typedef map< unsigned int, map<unsigned int, Edge> >::iterator vertex_it;
        typedef map<unsigned int, Edge>::iterator edge_it;

        // Iterate through all vertices
        for (vertex_it v_it = net.begin(); v_it != net.end(); v_it++) {
            map<unsigned int, Edge>& edges = v_it->second;
            
            // Iterate through all edges of a given vertex
            for (edge_it e_it = edges.begin(); e_it != edges.end(); e_it++) {
                Edge& edge = e_it->second;

                edge.flow = 0;
            }
        }
    }

    void add_flow(unsigned int orig, unsigned int dest, int flow) {
        check_edge(orig, dest);

        // ===== Check if flow value is at most the capacity =====
        stringstream ss;
        if ( (net[orig][dest].flow + flow) > net[orig][dest].capacity ) {
            ss << "Flow value " << (net[orig][dest].flow + flow) << 
                  " cannot be greater than capacity " << 
                  net[orig][dest].capacity << " on edge (" << orig << 
                  "," << dest << ")";
            throw logic_error( ss.str() );
        }

        if ( (net[dest][orig].flow - flow) > net[dest][orig].capacity ) {
            ss << "Flow value " << (net[dest][orig].flow - flow) << 
                  " cannot be greater than capacity " << 
                  net[orig][dest].capacity << " on edge (" << orig <<
                  "," << dest << ")";
            throw logic_error( ss.str() );
        }
        // =======================================================

        net[orig][dest].flow += flow;
        net[dest][orig].flow -= flow;
    }

    // ========================================================================



    // ========================== GETTERS & SETTERS ===========================
    
    int get_flow(unsigned int orig, unsigned int dest) {
        check_edge(orig, dest);    
        return net[orig][dest].flow;
    }


    int get_capacity(unsigned int orig, unsigned int dest) {
        return net[orig][dest].capacity;
    }


    map< unsigned int, map<unsigned int, Edge> > get_net() {
        return net;
    }

    // ========================================================================



    // =========================== GRAPH MODIFIERS ============================
    
    void add_vertex(unsigned int id) {
        if (net.find(id) != net.end()) {
            stringstream ss;
            ss << "Vertex " << id << " already exists";
            throw logic_error( ss.str() );
        }

        net[id] = map<unsigned int, Edge>();
    }


    void add_edge(unsigned int orig, unsigned int dest, int capacity) {
        stringstream ss;
        // Check if both vertices exist
        if (net.find(orig) == net.end() || net.find(dest) == net.end()) {
            ss << "Edge (" << orig << "," << dest << 
                  ") must be between existing vertices";
            throw logic_error( ss.str() );
        }

        // Check if capacity is positive
        if (capacity < 0) {
            ss << "Capacity value " << capacity << " must be positive";
            throw domain_error( ss.str() );
        }

        // Get (and create, if needed) both forward and backward edges 
        // (orig,dest) and (dest,orig), respectively.
        Edge& forwards = net[orig][dest];
        net[dest][orig];

        // Increase forward edge capacity
        forwards.capacity += capacity;
    }

    // ========================================================================



    // ======================= MAX FLOW (EDMONDS-KARP) ========================
    
    // TODO: function
    void max_flow () {
        clean_flow();
        

    }

    // ========================================================================
};



// =============================== IO FUNCTIONS ===============================

istream& read_problem_instance(istream& stream, Network& network) {
    unsigned int num_vertices, num_edges;

    if (!stream)
        return stream;

    if (!(stream >> num_edges) || !(stream >> num_vertices))
        throw runtime_error("Error reading number of edges and vertices");

    // Add vertices to the network
    for (unsigned int i = 1; i <= num_vertices; i++)
        network.add_vertex(i);

    // Add edges to the network
    for (unsigned int i = 0; i < num_edges; i++) {
        unsigned int origin, destination, capacity;
        
        // Read the line's data into the variables
        if (!(stream >> origin) || !(stream >> destination) || 
              !(stream >> capacity)) {
            stringstream ss;
            ss << "Error reading edge #" << (i + 1);
            throw runtime_error( ss.str() );
        }

        network.add_edge(origin, destination, capacity);
    }

    return stream;
}


void solve_multiple_instances() {
    fstream input_file;
    input_file.open("in.txt", fstream::in);

    if (!input_file)
        throw runtime_error("Could not open input file");

    unsigned int num_instances;
    if ( !(input_file >> num_instances) )
        throw runtime_error("Error reading number of instances");

    for (unsigned int i = 0; i < num_instances; i++) {
        try {
            Network network;        
            read_problem_instance(input_file, network);
        
            map< unsigned int, map<unsigned int, Edge> > net = network.get_net();

            cout << "Edges: " << endl;
            for (map< unsigned int, map<unsigned int, Edge> >::const_iterator it = net.begin(); it != net.end(); it++) {
                unsigned int id = it->first;
                map<unsigned int, Edge> edges = it->second;

                for (map<unsigned int, Edge>::const_iterator it2 = edges.begin(); it2 != edges.end(); it2++) {
                    cout << '(' << id << ',' << it2->first << ") => " << it2->second.flow << '/' << it2->second.capacity << endl;
                }
            }

            // TODO: set source and sink
        
        } catch (exception& e) {
            cout << e.what() << endl;
        }
    }
}

// ============================================================================



int main() {
    solve_multiple_instances();
    
    return 0;
}
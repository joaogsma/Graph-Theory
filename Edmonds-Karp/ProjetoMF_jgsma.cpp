#include <algorithm>
#include <exception>
#include <fstream>
#include <iostream>
#include <climits>
#include <unordered_map>
#include <queue>
#include <set>
#include <stack>
#include <stdexcept>
#include <string>
#include <sstream>
#include <utility>
#include <vector>

using std::cout;            using std::unordered_map;
using std::domain_error;    using std::pair;
using std::endl;            using std::queue;
using std::exception;       using std::runtime_error;
using std::fstream;         using std::set;
using std::getline;         using std::stack;
using std::istream;         using std::stringstream;
using std::logic_error;     using std::string;
using std::min;             using std::vector;

struct Max_Flow_Info {
    int max_flow;
    unordered_map< unsigned int, unordered_map<unsigned int, int> > min_cut_edges;
};

struct Edge {
    int flow, capacity;

    Edge() : flow(0), capacity(0) {}
    Edge(int flow, int capacity) : flow(flow), capacity(capacity) {}
};

class Network {
private:
    unordered_map< unsigned int, unordered_map<unsigned int, Edge> > net;
    unsigned int source, sink;

    // ======================= MISCELLANEOUS FUNCTIONS ========================

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

    // ========================================================================



    // ==================== MAX FLOW (EDMONDS-KARP) - BFS =====================
    
    int bfs_find_residual_path(vector<unsigned int>& path, 
          vector< pair<unsigned int, unsigned int> >& minimum_cut) {
        // ========== Error checks ==========
        // Empty graph
        if (net.size() == 0)
            return false;

        // Source and sink not initialized
        if (source == sink)
            return false;
        // ====================================


        // ========== Setup ==========
        // Erase possible values in variables
        path.clear();
        minimum_cut.clear();

        // Set up variables
        queue<unsigned int> vertex_order;
        unordered_map<unsigned int, unsigned int> parent_map;
        vertex_order.push(source);
        vector< pair<unsigned int, unsigned int> > saturated_edges;
        // ===========================


        // ========== BFS loop ==========
        while ( !vertex_order.empty() ) {
            unsigned int vertex = vertex_order.front();

            // Check if the a path was found
            if (vertex == sink)
                break;

            unordered_map<unsigned int, Edge>& edges = net[vertex];

            for (unordered_map<unsigned int, Edge>::const_iterator it = edges.begin(); 
                  it != edges.end(); it++) {
                unsigned int child = it->first;
                const Edge& edge = it->second;

                // True if the child has not been visited yet
                bool valid_child = parent_map.find(child) == parent_map.end();
                // True if there is available residual capacity
                bool valid_edge = (edge.capacity - edge.flow) > 0;
                // True if it is a saturated edge
                bool saturated_edge = edge.capacity > 0 && 
                      edge.capacity == edge.flow;

                // If this is a valid edge to a valid child, add the child to
                // the queue and store vertex as its parent
                if ( valid_child && valid_edge ) {
                    parent_map[child] = vertex;
                    vertex_order.push(child);
                } else if (valid_child && saturated_edge) {
                    saturated_edges.push_back( 
                          pair<unsigned int, unsigned int>(vertex, child) );
                }
            }

            vertex_order.pop();
        }
        // ==============================


        // ========== Path processing ==========
        // Check if a path was found
        if ( parent_map.find(sink) != parent_map.end() ) {
            stack<unsigned int> stack;
            stack.push(sink);

            // Add path vertices to stack in order to invert their order
            while (stack.top() != source)
                stack.push( parent_map[stack.top()] );

            int min_flow = INT_MAX;

            // Add stack vertices to path vector, in the correct order
            while ( !stack.empty() ) {
                unsigned int current = stack.top();
                stack.pop();
                path.push_back(current);

                // Keep track of the smallest residual capacity
                if (!stack.empty()) {
                    Edge& edge = net[current][stack.top()];
                    int residual_capacity = edge.capacity - edge.flow;
                    min_flow = min(min_flow, residual_capacity);
                }
            }

            return min_flow;
        } else {
            minimum_cut.insert(minimum_cut.end(), saturated_edges.begin(), 
                  saturated_edges.end());
        }
        // =====================================


        return 0;
    }

    // ========================================================================



public:
    // ====================== FLOW MODIFYING FUNCTIONS ========================

    void clean_flow() {
        typedef unordered_map< unsigned int, unordered_map<unsigned int, Edge> >::iterator vertex_it;
        typedef unordered_map<unsigned int, Edge>::iterator edge_it;

        // Iterate through all vertices
        for (vertex_it v_it = net.begin(); v_it != net.end(); v_it++) {
            unordered_map<unsigned int, Edge>& edges = v_it->second;
            
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


    void set_source_and_sink(unsigned int source, unsigned int sink) {
        if (source == sink)
            throw logic_error("Source and Sink must be different vertices");

        stringstream ss;
        if (net.find(source) == net.end()) {
            ss << "Vertex " << source << " does not exist";
            throw logic_error( ss.str() );
        }

        if (net.find(sink) == net.end()) {
            ss << "Vertex " << sink << " does not exist";
            throw logic_error( ss.str() );
        }

        this->source = source;
        this->sink = sink;
    }

    // ========================================================================



    // ======================= GRAPH BUILDING FUNCTIONS =======================
    
    void add_vertex(unsigned int id) {
        if (net.find(id) != net.end()) {
            stringstream ss;
            ss << "Vertex " << id << " already exists";
            throw logic_error( ss.str() );
        }

        net[id] = unordered_map<unsigned int, Edge>();
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
    
    void max_flow (Max_Flow_Info& info) {
        clean_flow();
        
        int total_flow = 0;
        vector<unsigned int> path;
        vector< pair<unsigned int, unsigned int> > minimum_cut;
        int min_flow;

        // Iterate while there is a path between source and skin
        min_flow = bfs_find_residual_path(path, minimum_cut);
        while ( min_flow ) {
            total_flow += min_flow;

            // Add min_flow to all edges in the returned path
            for (vector<unsigned int>::size_type i = 0; i < path.size() - 1; i++)
                add_flow( path[i], path[i+1], min_flow);

            // Try to find another path
            min_flow = bfs_find_residual_path(path, minimum_cut);
        }

        // Store maximum flow value
        info.max_flow = total_flow;

        // Store the saturated edges
        typedef vector< pair<unsigned int, unsigned int> > edge_vec;
        for (edge_vec::const_iterator it = minimum_cut.begin(); 
              it != minimum_cut.end(); it++) {
            unsigned int orig = it->first, dest = it->second;
            info.min_cut_edges[orig][dest] = net[orig][dest].flow;
        }
    }

    // ========================================================================
};



// =============================== IO FUNCTIONS ===============================

string instance_to_string(int instance, const Max_Flow_Info& info) {
    stringstream ss;

    ss << "O fluxo máximo na rede " << instance << " tem valor " << 
          info.max_flow << ". Um corte saturado na rede " << instance << 
          " com este fluxo contém as arestas: " << endl;

    for (unordered_map< unsigned int, unordered_map<unsigned int, int> >::const_iterator orig_it = info.min_cut_edges.begin(); 
          orig_it != info.min_cut_edges.end(); orig_it++) {
        unsigned int orig = orig_it->first;
        const unordered_map<unsigned int, int>& edges = orig_it->second;
            
        for (unordered_map<unsigned int, int>::const_iterator dest_it = edges.begin();
              dest_it != edges.end(); dest_it++) {
            unsigned int dest = dest_it->first;
            int flow = dest_it->second;
            ss << orig << "," << dest << " " << flow << endl;
        }
    }

    return ss.str();
}


istream& read_problem_instance(istream& stream, Network& network, int instance) {
    unsigned int num_vertices, num_edges;

    if (!stream)
        return stream;

    if (!(stream >> num_edges) || !(stream >> num_vertices)) {
        stringstream ss;
        ss << "On instance " << instance << 
              " => Could not read number of edges and vertices";
        throw runtime_error( ss.str() );
    }

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
            ss << "On instance " << instance << " => Could not read edge #" << (i + 1);
            throw runtime_error( ss.str() );
        }

        network.add_edge(origin, destination, capacity);
    }

    // Set network source and sink
    if (num_vertices > 0) 
        network.set_source_and_sink(1, num_vertices);

    return stream;
}


void solve_multiple_instances(string in_filename, string out_filename) {
    fstream input_file, output_file;
    input_file.open(in_filename, fstream::in);
    output_file.open(out_filename, fstream::out);

    if (!input_file)
        throw runtime_error("Could not open file " + in_filename);

    if (!output_file)
        throw runtime_error("Could not open file " + out_filename);        

    unsigned int num_instances;
    if ( !(input_file >> num_instances) )
        throw runtime_error("Could not read number of instances");

    for (unsigned int i = 0; i < num_instances; i++) {
        Network network;
        read_problem_instance(input_file, network, i+1);

        Max_Flow_Info info;
        network.max_flow(info);

        string output = instance_to_string(i+1, info);

        if (i != 0)
            if ( !(output_file << endl) )
                throw runtime_error("Could not write to file " + out_filename);        

        if ( !(output_file << output) )
            throw runtime_error("Could not write to file " + out_filename);
    }

    input_file.close();
    output_file.close();
}

// ============================================================================



int main() {
    try {
        solve_multiple_instances("ProjetoMF_input.txt", 
              "ProjetoMF_output_jgsma.txt");
    } catch (exception& e) {
        cout << "ERROR: " << e.what() << endl;
    }
    
    return 0;
}
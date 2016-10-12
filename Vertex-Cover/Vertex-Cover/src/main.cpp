#include <algorithm>
#include <chrono>
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

using std::chrono::duration;		using std::mutex;
using std::chrono::duration_cast;	using std::pair;
using std::chrono::steady_clock;	using std::random_device;
using std::cout;					using std::runtime_error;
using std::endl;					using std::sort;
using std::fstream;					using std::string;
using std::ostream;					using std::stringstream;
using std::logic_error;				using std::set;
using std::map;						using std::thread;
using std::mt19937;					using std::vector;


typedef unsigned int vertex_id;

static unsigned int finished_cout = 0;
static const unsigned int num_repetitions = 300;
static fstream output_file;
static random_device rd;
static mt19937 mt( rd() );
static mutex lock;

// ============================================================================
// ============================ CLASS DEFINITIONS =============================
// ============================================================================

class Simple_Graph {
public:
	const map<vertex_id, set<vertex_id> >& edges() const { return graph; }
    
	void clear() { graph.clear(); }

    void clear_edges();
    
    void add_vertex(vertex_id id);

    void add_edge(vertex_id a, vertex_id b);

    void greedy_vertex_cover(set<vertex_id> &cover) const;

    void optimal_vertex_cover(set<vertex_id> &cover) const;

    string to_string() const;

	string to_string_output() const;

    void generate_as_random(unsigned int num_vertices, unsigned int num_edges);

private:
    map<vertex_id, set<vertex_id> > graph;

    bool is_valid_cover(const set<vertex_id> &combination) const;
    
    vertex_id greatest_degree_vertex(
        const map<vertex_id, set<vertex_id> > &graph) const;

    vertex_id pick_random(const set<vertex_id>& vertex_set) const;
    
	void add_random_edges(const set<vertex_id>& vertex_set, 
        unsigned int num_new_edges);
    
	vector<vertex_id> loop_erasing_random_walk(
        const set<vertex_id>& tree_vertices, const set<vertex_id>& vertex_set) const;

    void wilson_spanning_tree(const set<vertex_id>& vertex_set);

};


struct Analysis_Info {
	Simple_Graph graph;
	set<vertex_id> greedy_cover;
	set<vertex_id> optimal_cover;
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


vertex_id Simple_Graph::greatest_degree_vertex(
    const map<vertex_id, set<vertex_id> > &graph) const
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
	/*	Store all the vertices that can still support new edges. A vertex will be 
		removed from this set if it already has |V|-1 edges, i.e., it has an edge 
		to every other vertex */
	set<vertex_id> available_vertices = vertex_set;

    while( num_new_edges > 0 )
    {

        // ===== Pick two different vertices =====
        vertex_id vertex1 = pick_random(available_vertices);
        
        /*	Compute the set of available vertices not adjacent to vertex1. 
			These are the options that can have an edge to vertex1 */
        set<vertex_id> not_adjacent = available_vertices;
        not_adjacent.erase( vertex1 );
        set<vertex_id>& vertex1_adjacent = graph[vertex1];
		
		for (set<vertex_id>::const_iterator it = vertex1_adjacent.begin();
			it != vertex1_adjacent.end(); it++)
		{
			// Check if the vertex *it is an available one
			set<vertex_id>::iterator pos = not_adjacent.find( *it );
			/*	If the vertex *it is no longer available, remove it from 
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
	unsigned int num_edges)
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

    // Add the specified number of vertices to both the graph and the vertex set
    for (unsigned int i = 1; i <= num_vertices; i++)
    {
        add_vertex(i);
        vertex_set.insert(i);
    }

    // Generate a random spanning tree to ensure the graph is connected
    wilson_spanning_tree(vertex_set);

    /*	Add random edges to achieve the specified density. So far, there are
		num_vertices-1 edges, since the graph is a tree */
    add_random_edges(vertex_set, num_edges - (num_vertices-1));
}

// ============================================================================



// ============================================================================
// ========================= VERTEX COVER FUNCTIONS =========================
// ============================================================================

void Simple_Graph::greedy_vertex_cover(set<vertex_id> &cover) const
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


void Simple_Graph::optimal_vertex_cover(set<vertex_id> &cover) const
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

bool compare(const Analysis_Info& a, const Analysis_Info& b)
{
	double opt_dist_a = static_cast<double>(a.greedy_cover.size()) /
		a.optimal_cover.size();
	double opt_dist_b = static_cast<double>(b.greedy_cover.size()) /
		b.optimal_cover.size();

	return opt_dist_a > opt_dist_b;
}


struct Thread_Params {
	unsigned int num_vertices, num_edges;
	vector<Analysis_Info> *worst_graphs;
	double *sum;
};

void run_iteration(Thread_Params params) {
	// Initialize the graph randomly
	Simple_Graph g;
	g.generate_as_random(params.num_vertices, params.num_edges);

	// Compute the optimal vertex cover through brute-force
	set<vertex_id> optimal_cover;
	g.optimal_vertex_cover(optimal_cover);

	// Compute the vertex cover though the greedy approach
	set<vertex_id> greedy_cover;
	g.greedy_vertex_cover(greedy_cover);

	double opt_dist = static_cast<double>(greedy_cover.size()) / optimal_cover.size();

	Analysis_Info analysis;
	analysis.graph = g;
	analysis.greedy_cover = greedy_cover;
	analysis.optimal_cover = optimal_cover;

	// ===== Exclusive Zone =====
	lock.lock();
	
	cout << ++finished_cout << endl;

	(*params.sum) += opt_dist;
	params.worst_graphs->push_back(analysis);

	// Keep only the 3 worst seen so far
	if (params.worst_graphs->size() > 3)
	{
		sort(params.worst_graphs->begin(), params.worst_graphs->end(), compare);
		params.worst_graphs->resize(3);
	}
	
	lock.unlock();
	// ==========================
}


double analyze_category(unsigned int num_vertices,
	unsigned int num_edges, vector<Analysis_Info>& worst_graphs)
{
	worst_graphs.clear();

	double sum = 0;

	vector<thread> thread_pool;

	for (unsigned int iteration = 0; iteration < num_repetitions; iteration++)
	{
		Thread_Params params;
		params.num_vertices = num_vertices;
		params.num_edges = num_edges;
		params.sum = &sum;
		params.worst_graphs = &worst_graphs;
		
		thread_pool.push_back( thread(run_iteration, params) );

		if (thread_pool.size() == 8)
		{
			// Wait for the 8 current threads
			for (int i = 0; i < thread_pool.size(); i++)
				thread_pool[i].join();

			thread_pool.clear();
		}
	}

	// Wait for remaining threads
	for (int i = 0; i < thread_pool.size(); i++)
		thread_pool[i].join();

	finished_cout = 0;

	return sum / num_repetitions;
}


void run_analysis(const vector<pair<unsigned int, unsigned int> >& configurations)
{
	typedef vector<pair<unsigned int, unsigned int>>::const_iterator config_c_it;
	for (config_c_it it = configurations.begin(); it != configurations.end(); it++)
	{
		unsigned int num_vertices = it->first;
		unsigned int num_edges = it->second;
		
		vector<Analysis_Info> worst_graphs;
		double average = analyze_category(num_vertices, num_edges, worst_graphs);

		double worst_opt_dist = static_cast<double>(worst_graphs[0].greedy_cover.size()) /
			worst_graphs[0].optimal_cover.size();
		
		stringstream ss;

		// ===== Print configuration and worst distance to optimal =====
		ss << "========== Configuração: |V| = " << num_vertices << ", |E| = " << 
			num_edges << " ==========" << endl;
		ss << "Distância média para o ótimo: " << average << endl;
		ss << "Maior distância para o ótimo: " << worst_opt_dist << endl;
		// =============================================================

		ss << endl;

		for (vector<Analysis_Info>::const_iterator worst_it = worst_graphs.begin();
			worst_it != worst_graphs.end(); worst_it++)
		{
			ss << "Número de vértices = " << worst_it->graph.edges().size() << endl <<
				worst_it->graph.to_string_output() << endl;

			// ===== Print the greedy solution set to the output stream =====
			ss << "Solução Gulosa: {";
			for (set<vertex_id>::const_iterator greedy_sol_it = worst_it->greedy_cover.begin();
				greedy_sol_it != worst_it->greedy_cover.end(); greedy_sol_it++)
			{
				if (greedy_sol_it != worst_it->greedy_cover.begin())
					ss << ", ";
				ss << *greedy_sol_it;
			}
			ss << "}" << endl;
			// ==============================================================


			// ===== Print the optimal solution set to the output stream =====
			ss << "Solução Ótima: {";
			for (set<vertex_id>::const_iterator opt_sol_it = worst_it->optimal_cover.begin();
				opt_sol_it != worst_it->optimal_cover.end(); opt_sol_it++)
			{
				if (opt_sol_it != worst_it->optimal_cover.begin())
					ss << ", ";
				ss << *opt_sol_it;
			}
			ss << "}" << endl;
			// ===============================================================

			ss << endl;

		}

		// Write to file, throw runtime_error if cannot write
		if (!(output_file << ss.str()))
			throw runtime_error("Could not write to output file");
		
		cout << "Finished a category" << endl;

		output_file.flush();
	}
}

// ============================================================================



// ============================================================================
// ============================== MAIN FUNCTION ===============================
// ============================================================================


int main()
{
	try {
		vector< pair<unsigned int, unsigned int> > configurations;
		configurations.push_back(pair<unsigned int, unsigned int>(15, 30));
		configurations.push_back(pair<unsigned int, unsigned int>(15, 45));
		configurations.push_back(pair<unsigned int, unsigned int>(15, 45));
		configurations.push_back(pair<unsigned int, unsigned int>(15, 57));
		
		configurations.push_back(pair<unsigned int, unsigned int>(20, 40));
		configurations.push_back(pair<unsigned int, unsigned int>(20, 60));
		configurations.push_back(pair<unsigned int, unsigned int>(20, 80));
		configurations.push_back(pair<unsigned int, unsigned int>(20, 100));
		
		
		output_file.open("GVC_jgsma_output.txt", fstream::out);
		
		steady_clock::time_point begin_time = steady_clock::now();

		run_analysis(configurations);

		steady_clock::time_point end_time = steady_clock::now();

		int elapsed_time_hours = duration_cast<std::chrono::hours>(
			end_time - begin_time).count();

		int elapsed_time_minutes = duration_cast<std::chrono::minutes>(
			end_time - begin_time).count() - elapsed_time_hours*60;

		int elapsed_time_seconds = duration_cast<std::chrono::seconds>(
			end_time - begin_time).count() - elapsed_time_hours*3600 - 
			elapsed_time_minutes*60;

		cout << "TEMPO DE EXECUÇÃO: " << elapsed_time_hours << " horas, " << 
			elapsed_time_minutes << " minutos e "  << elapsed_time_seconds << 
			" segundos" << endl;

		output_file.close();

		system("pause");
	}
	catch (std::exception& e)
	{
		cout << "ERROR:  " << e.what() << endl;
	}

	return 0;
}
// ============================================================================
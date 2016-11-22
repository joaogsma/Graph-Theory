#ifndef __GUARD_MAIN__
#define __GUARD_MAIN__

#include <algorithm>
#include <climits>
#include <chrono>
#include <fstream>
#include <iterator>
#include <map>
#include <numeric>
#include <vector>
#include <random>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>

#include "gurobi_c++.h"

using std::accumulate;
using std::back_inserter;
using std::chrono::steady_clock;
using std::chrono::duration_cast;
using std::copy;
using std::cout;
using std::endl;
using std::ifstream;
using std::istream;
using std::make_pair;
using std::map;
using std::min_element;
using std::mt19937;
using std::ofstream;
using std::ostream;
using std::pair;
using std::random_device;
using std::runtime_error;
using std::set;
using std::string;
using std::stringstream;
using std::transform;
using std::uniform_int_distribution;
using std::vector;

// ============================================================================
// ============================== CONFIGURATIONS ==============================
// ============================================================================

static const int min_distance = 1;
static const int max_distance = 10;

static const int min_client_demand = 10;
static const int max_client_demand = 50;

static const int min_opening_cost = 80;
static const int max_opening_cost = 300;

static const int min_clients_per_facility = 5;
static const int max_clients_per_facility = 10;

static const int min_facilities_per_client = 3;
static const int max_facilities_per_client = 4;

static const int artificial_inf = std::max(1000 * max_distance, 99999);

// ============================================================================


static random_device rd;
static mt19937 mt_engine( rd() );
static GRBEnv env;


// ============================================================================
// ============================== LINEAR PROGRAM ==============================
// ============================================================================

// ============= CLASS DECLARATION =============
/*  The Linear_Program class assumes a linear program in standard form, i.e., 
    the objective function is to be maximized, all inequalities are in the 
    form ... <= ... and all variables in th e solution must be >= 0. */
class Linear_Program {
friend int main();
public:
    Linear_Program(const vector<double>& objective_vector, 
        const vector<vector<double> >& constraints,
        const vector<char>& constraint_types): objective_vector(objective_vector), 
        constraints(constraints), constraint_types(constraint_types) { validate(); }

    void solve(vector<double>& solution, double& value, bool do_maximize = true, 
        bool relaxed = false) const;

private:
    vector<double> objective_vector;    // Coefficients of the objective function to be maximized
    vector<vector<double> > constraints;   // A matrix of the linear program in standard form
    vector<char> constraint_types;            // b vector of the linear program in standard form

    void validate() const;

    void gurobi_solve(vector<double>& solution, double& max, bool do_maximize, 
        bool relaxed) const;
};
// =============================================


// ============= CLASS DEFINITION =============
void Linear_Program::validate() const
{
    const string& error_message = "Invalid Linear Programming configuration";

    // Check if either vector is empty
    if (objective_vector.empty() || constraints.empty() || constraint_types.empty())
        throw runtime_error(error_message);
    
    // Check if objective vector of coefficients contains coefficients 
    // for all variables
    if (objective_vector.size() != constraints[0].size() - 1)
        throw runtime_error(error_message);

    // Check if there is one type value for each equation
    if (constraint_types.size() != constraints.size())
        throw runtime_error(error_message);

    // Check if all constraints contain coefficients for all variables
    vector<double>::size_type eq_size = constraints[0].size();
    for (vector<double>::size_type i = 1; i < constraints.size(); ++i)
        if (constraints[i].size() != eq_size) throw runtime_error(error_message);
}

void Linear_Program::solve(vector<double>& solution, double& value, 
    bool do_maximize, bool relaxed) const
{
    gurobi_solve(solution, value, do_maximize, relaxed);
}

void Linear_Program::gurobi_solve(vector<double>& solution, double& value, 
    bool do_maximize, bool relaxed) const
{
    size_t num_variables = objective_vector.size();

    // Create the gurobi model
    GRBModel model(env);

    // Add the variables to the model
    vector<GRBVar> variables;
    for (size_t i = 0; i < num_variables; ++i)
    {
        variables.push_back(model.addVar(0., (relaxed ? GRB_INFINITY : 1.),
            0., GRB_CONTINUOUS));
    }

    // Add the objective function to the model
    GRBLinExpr obj_expr;
    obj_expr.addTerms(&objective_vector[0], &variables[0], num_variables);
    model.setObjective(obj_expr, (do_maximize ? GRB_MAXIMIZE : GRB_MINIMIZE));

    // Add the constraints to the model
    for (size_t i = 0; i < constraints.size(); ++i)
    {
        GRBLinExpr expr;
        expr.addTerms(&constraints[i][0], &variables[0], num_variables);
        model.addConstr(expr, (constraint_types[i] == '=' ? GRB_EQUAL : GRB_LESS_EQUAL),
            constraints[i].back());
    }

    model.optimize();

    solution.clear();
    for (size_t i = 0; i < num_variables; ++i)
        solution.push_back(variables[i].get(GRB_DoubleAttr_X));

    value = model.get(GRB_DoubleAttr_ObjVal);
}

// ============================================

// ============================================================================



// ============================================================================
// ================================== GRAPH ===================================
// ============================================================================

struct Graph {
    vector<vector<int> > distance_cost;
    vector<int> facility_costs;
    vector<int> client_demands;

    int facilities() const { return (int)facility_costs.size(); }

    int clients() const { return (int)client_demands.size(); }

private:
    static int get_key(const pair<int, int>& map_it) { return map_it.first; }
};

Graph random_graph(int num_facilities, int num_clients)
{
    Graph problem;

    // Initially populate the distances with 0 distances
    problem.distance_cost = vector<vector<int> >(num_facilities,
        vector<int>(num_clients, 0));

    uniform_int_distribution<int> fac_cost_dist(min_opening_cost, max_opening_cost);

    // Initialize the opening costs with random values from [min_opening_cost, max_opening_cost]
    for (int facility = 0; facility < num_facilities; ++facility)
        problem.facility_costs.push_back(fac_cost_dist(mt_engine));

    uniform_int_distribution<int> client_dem_dist(min_client_demand, max_client_demand);

    // Initialize the client demands with random values from [min_client_demand, max_client_demand]
    for (int client = 0; client < num_clients; ++client)
        problem.client_demands.push_back(client_dem_dist(mt_engine));

    // Assign the minimum number of facilities for each client
    for (int client = 0; client < num_clients; ++client)
    {
        int has_facilities = 0;
        vector<int> available;
        for (int i = 0; i < num_facilities; ++i) available.push_back(i);

        while (has_facilities < min_facilities_per_client)
        {
            if (available.empty())    // Failed, try again
                return random_graph(num_facilities, num_clients);

            uniform_int_distribution<int> facility_dist(0, available.size() - 1);
            int facility_idx = facility_dist(mt_engine);
            int facility = available[facility_idx];    // Random facility

                                                       // Remove facility from the available vector
            available.erase(available.begin() + facility_idx);

            /*  Mark this client as supported by the facility if it doest not
            already support the maximum number of clients */
            if (accumulate(problem.distance_cost[facility].begin(),
                problem.distance_cost[facility].end(), 0) < max_clients_per_facility)
            {
                problem.distance_cost[facility][client] = 1;
                has_facilities++;
            }
        }
    }

    // Assign the (remaining) minimum number of clients to each facility
    for (int facility = 0; facility < num_facilities; ++facility)
    {
        int has_clients = accumulate(problem.distance_cost[facility].begin(),
            problem.distance_cost[facility].end(), 0);

        vector<int> available;
        for (int i = 0; i < num_clients; ++i) available.push_back(i);

        while (has_clients < min_clients_per_facility)
        {
            if (available.empty())    // Failed, try again
                return random_graph(num_facilities, num_clients);

            uniform_int_distribution<int> client_dist(0, available.size() - 1);
            int client_idx = client_dist(mt_engine);
            int client = available[client_idx];

            // Remove client from the available vector
            available.erase(available.begin() + client_idx);

            int client_facilities = 0;
            for (int i = 0; i < num_facilities; ++i)
                client_facilities += problem.distance_cost[i][client];

            /*  Support client with this facility if the client is not yet
            supported by the max number of facilities */
            if (client_facilities < max_facilities_per_client)
            {
                problem.distance_cost[facility][client] = 1;
                has_clients++;
            }
        }
    }

    int possible_extra_conections = (max_facilities_per_client - min_facilities_per_client) *
        (max_clients_per_facility - min_clients_per_facility);

    uniform_int_distribution<int> extra_dist(0, possible_extra_conections);
    int extra_connections = extra_dist(mt_engine);

    // Vector containing all positions connections not already chosen
    vector<pair<int, int> > possible_connections;

    // Fill the possible_conections vector
    for (int facility = 0; facility < num_facilities; ++facility)
        for (int client = 0; client < num_clients; ++client)
            if (problem.distance_cost[facility][client] == 0)
                possible_connections.push_back(make_pair(facility, client));

    for (int i = 0; i < extra_connections && !possible_connections.empty(); ++i)
    {
        uniform_int_distribution<int> position_dist(0, possible_connections.size() - 1);
        int idx = position_dist(mt_engine);
        int facility = possible_connections[idx].first;
        int client = possible_connections[idx].second;

        possible_connections.erase(possible_connections.begin() + idx);

        // Count the number of facilities associated with this client
        int has_facilities = 0;
        for (int i = 0; i < num_facilities; ++i)
            has_facilities += problem.distance_cost[i][client];
        // Client cannot have more facilities
        if (has_facilities == max_facilities_per_client) continue;

        // Count the number of clients associated with this facility
        int has_clients = 0;
        for (int i = 0; i < num_clients; ++i)
            has_clients += problem.distance_cost[facility][i];
        // Facility cannot have more clients
        if (has_clients == max_clients_per_facility) continue;

        problem.distance_cost[facility][client] = 1;
    }


    uniform_int_distribution<int> distance_dist(min_distance, max_distance);

    // Substitute all 0s by "infinite" (no connection) and all 1s by a random distance
    for (int facility = 0; facility < num_facilities; ++facility)
    {
        for (int client = 0; client < num_clients; ++client)
        {
            int& value = problem.distance_cost[facility][client];
            value = value ? distance_dist(mt_engine) : artificial_inf;
        }
    }

    return problem;
}

// ============================================================================



// ============================================================================
// ================================= UTILITY ==================================
// ============================================================================

inline int edge_pos_primal(int client, int facility, int num_clients, 
    int num_facilities = 0)
{
    return num_facilities + (facility * num_clients) + client;
}

inline int edge_pos_dual(int facility, int client, int num_facilities, 
    int num_clients = 0)
{
    return num_clients + (client * num_facilities) + facility;
}

// ============================================================================



// ============================================================================
// =============================== BRUTE-FORCE ================================
// ============================================================================

inline bool is_open(unsigned long long configuration, int facility)
{
    return ((1ll << facility) & configuration) > 0;
}

void optimal(const Graph& problem, map<int, int>& client_assignments, 
    map<int, vector<int> >& facility_assignments, double& min_cost)
{
    client_assignments.clear();
    facility_assignments.clear();

    int num_facilities = problem.facilities();
    int num_clients = problem.clients();
    int num_variables = num_facilities * num_clients;

    unsigned long long best_opened;
    vector<double> best_solution;
    double best_cost = INT_MAX;

    unsigned long long out_of_range = 1ll << num_facilities;

    /*  Test all combinations of facilities. Each combination is codified in an
    unsigned long long variable, where a 1 value in the nth bit signals that
    the facility in the nth position of the facilities vector is active */
    for (unsigned long long combination = 1; combination < out_of_range; ++combination)
    {
        /*  Objecive function to be minimized. The coefficients will all be negated
            before they are added, since a LP in standard form requires a function
            to be maximized */
        vector<double> objective_vector(num_variables, 0);

        // Copy facility-client edges (distance cost * client demand)
        for (int facility = 0; facility < num_facilities; ++facility)
        {
            for (int client = 0; client < num_clients; ++client)
            {
                int pos = edge_pos_primal(client, facility, num_clients);
                objective_vector[pos] = problem.distance_cost.at(facility).at(client) *
                    problem.client_demands[client];
            }
        }

        vector<vector<double> > constraints;
        vector<char> constraint_types;

        /*  Constraints guaranteeing that the demands of every client is satified.
            These specify that the summation of x_ij, through all facilities i,
            is equal to one, for all clients j */
        for (int client = 0; client < num_clients; ++client)
        {
            vector<double> constraint(num_variables+1, 0);
            
            for (int facility = 0; facility < num_facilities; ++facility)
            {
                int position = edge_pos_primal(client, facility, num_clients);
                constraint[position] = 1;
            }
            constraint.back() = 1;

            constraints.push_back(constraint);
            constraint_types.push_back('=');
        }

        /*  Constraints guaranteeing that clients are supplied only from open
            facilities. These specify that x_ij <= y_j for all facility i and
            client j */
        for (int facility = 0; facility < num_facilities; ++facility)
        {
            for (int client = 0; client < num_clients; ++client)
            {
                vector<double> constraint(num_variables+1, 0);

                int pos = edge_pos_primal(client, facility, num_clients);
                constraint[pos] = 1;
                constraint.back() = is_open(combination, facility);

                constraints.push_back(constraint);
                constraint_types.push_back('<');
            }
        }

        // Create and solve a LP problem in order to find the remaining variables
        Linear_Program lp(objective_vector, constraints, constraint_types);

        vector<double> solution;
        double min;
        double opening_cost = 0;
        
        lp.solve(solution, min, false);

        for (int i = 0; i < num_facilities; ++i)
            opening_cost += is_open(combination, i) * problem.facility_costs[i];

        double cost = min + opening_cost;

        if (cost < best_cost)
        {
            best_solution = solution;
            best_opened = combination;
            best_cost = cost;
        }
    }

    min_cost = best_cost;
    
    for (int facility = 0; facility < num_facilities; ++facility)
    {
        for (int client = 0; client < num_clients; ++client)
        {
            int position = edge_pos_primal(client, facility, num_clients);
            if ( best_solution[position] == 1 )
            {
                client_assignments[client] = facility;
                facility_assignments[facility].push_back(client);
            }
        }
    }

}

// ============================================================================



// ============================================================================
// ================================ HEURISTICS ================================           
// ============================================================================

struct Cluster {
    set<int> clients, facilities;
};

inline double objective_fn_contribution(const Graph& problem, int facility,
    const set<int>& clients)
{
    auto lambda_acc_client_cost = [&](double acc, int client) -> double 
    {
        int distance = problem.distance_cost.at(facility).at(client);
        if (distance == artificial_inf) distance = 0;
        return acc + distance * problem.client_demands[client]; 
    };
    
    return accumulate(clients.begin(), clients.end(), 
        (double) problem.facility_costs[facility], lambda_acc_client_cost);
}

inline int find_cluster_seed(const map<int, int>& client_assignments,
    const vector<double>& dual_solution, int num_clients)
{
    double min_val = INT_MAX;
    int client_seed;

    for (int client = 0; client < num_clients; ++client)
    {
        // Ignore clients that were already assigned to clusters
        if ( client_assignments.find(client) != client_assignments.end() )
            continue;

        /*  If this is a client with a smaller vj than the best one so far, 
            update the variables */
        if ( dual_solution[client] < min_val )
        {
            min_val = dual_solution[client];
            client_seed = client;
        }
    }

    // Signals that no unassigned client exists
    if (min_val == INT_MAX) return -1;

    return client_seed;
}

inline void formulate_primal(const Graph& problem, vector<double>& objective,
    vector<vector<double> >& constraints, vector<char>& constraint_types)
{
    int num_facilities = problem.facilities();
    int num_clients = problem.clients();
    int num_variables = num_facilities + num_facilities * num_clients;

    // ========== Set the objective function ==========
    objective.resize(num_variables);
    
    // Copy the facility costs (coefficients of the first num_facilities variables)
    copy(problem.facility_costs.begin(), problem.facility_costs.end(),
        objective.begin());

    // Set the coefficients of the xij variables
    for (int facility = 0; facility < num_facilities; ++facility)
    {
        for (int client = 0; client < num_clients; ++client)
        {
            int pos = edge_pos_primal(client, facility, num_clients, num_facilities);
            
            /*  Set the value of xij to be dij*wj, where dij is the cost per product 
                unit and wj is the client demand, in product units */
            objective[pos] = problem.distance_cost.at(facility).at(client) *
                problem.client_demands[client];
        }
    }
    // ================================================


    // ========== Set the sum(xij) = 1 constraints ==========
    for (int client = 0; client < num_clients; ++client)
    {
        vector<double> constraint(num_variables + 1, 0.);
        
        for (int facility = 0; facility < num_facilities; ++facility)
        {
            int pos = edge_pos_primal(client, facility, num_clients, num_facilities);
            // xij coefficient is 1
            constraint[pos] = 1.;
        }
        
        // Sum of xij, through all facilities i, must be equal to 1
        constraint.back() = 1.;

        constraints.push_back( constraint );
        constraint_types.push_back( '=' );
    }
    // ======================================================


    // ========== Set the xij <= yi constraints ==========
    for (int facility = 0; facility < num_facilities; ++facility)
    {
        for (int client = 0; client < num_clients; ++client)
        {
            vector<double> constraint(num_variables + 1, 0.);

            // yi position is set to -1
            constraint[facility] = -1.;
            // xij position is set to 1
            int pos = edge_pos_primal(client, facility, num_clients, num_facilities);
            constraint[pos] = 1.;

            constraints.push_back(constraint);
            constraint_types.push_back('<');
        }
    }
    // ===================================================
}

inline void formulate_dual(const Graph& problem, vector<double>& objective,
    vector<vector<double> >& constraints, vector<char>& constraint_types)
{
    int num_facilities = problem.facilities();
    int num_clients = problem.clients();
    int num_variables = num_clients + num_facilities * num_clients;

    // ========== Set the objective function ==========
    // Initial num_clients positions are set to 1
    objective = vector<double>(num_clients, 1.);
    /*  Expand the array to size num_variables. Positions 
        [num_clients, num_variables) are set to 0 */
    objective.resize(num_variables);
    // ================================================


    // ========== Set vj - Beta_ij <= cij constraints ==========
    for (int client = 0; client < num_clients; ++client)
    {
        for (int facility = 0; facility < num_facilities; ++facility)
        {
            vector<double> constraint(num_variables + 1, 0.);

            // Set vj variable to 1
            constraint[client] = 1.;
            // Set Beta_ij variable to -1
            int pos = edge_pos_dual(facility, client, num_facilities, num_clients);
            constraint[pos] = -1.;

            // Set the right-hand-side to dij*wij
            constraint.back() = problem.distance_cost.at(facility).at(client) * 
                problem.client_demands[client];

            constraints.push_back(constraint);
            constraint_types.push_back('<');
        }
    }
    // =========================================================


    // ========== Set sum(Beta_ij) <= fi constraints ==========
    for (int facility = 0; facility < num_facilities; ++facility)
    {
        vector<double> constraint(num_variables + 1, 0.);
        
        for (int client = 0; client < num_clients; ++client)
        {
            // Set position Beta_ij to 1
            int pos = edge_pos_dual(facility, client, num_facilities, num_clients);
            constraint[pos] = 1.;
        }

        // Set the right-hand-side to fi
        constraint.back() = problem.facility_costs[facility];

        constraints.push_back(constraint);
        constraint_types.push_back('<');
    }
    // ========================================================
}

void lp_rounding(const Graph& problem, map<int, int>& client_assignments, 
    map<int, vector<int> >& facility_assignments, double& min_cost, 
    bool modified = false)
{
    client_assignments.clear();
    facility_assignments.clear();

    // ========== Solve the primal and dual LP problems ==========
    vector<double> objective_primal, objective_dual;
    vector<vector<double> > constraints_primal, constraints_dual;
    vector<char> constraint_types_primal, constraint_types_dual;

    // Formulate LP problems
    formulate_primal(problem, objective_primal, constraints_primal, 
        constraint_types_primal);
    formulate_dual(problem, objective_dual, constraints_dual, 
        constraint_types_dual);

    // Create Linear_Program objects
    Linear_Program primal_lp(objective_primal, constraints_primal, constraint_types_primal);
    Linear_Program dual_lp(objective_dual, constraints_dual, constraint_types_dual);

    vector<double> solution_primal, solution_dual;
    double cost_primal, cost_dual;
    
    // Solve LP problems
    primal_lp.solve(solution_primal, cost_primal, false, true);
    dual_lp.solve(solution_dual, cost_dual, true, true);

    // ===========================================================


    // ========== Divide the clients into clusters ==========
    int num_facilities = problem.facilities();
    int num_clients = problem.clients();

    double cost = 0;

    int client_seed = find_cluster_seed(client_assignments, solution_dual, num_clients);

    // Create a new cluster and 
    while ( client_seed != -1 )
    {
        Cluster cluster;
        
        // Add to the cluster all facilities that serve the client seed
        for (int facility = 0; facility < num_facilities; ++facility)
        {
            // Add this facility to the cluster if the client seed is fractionally served by it
            int pos = edge_pos_primal(client_seed, facility, num_clients, num_facilities);
            if (solution_primal[pos] > 0) cluster.facilities.insert(facility);
        }

        // Add to the cluster all clients served by a facility in the cluster
        for (set<int>::const_iterator facility_it = cluster.facilities.begin();
            facility_it != cluster.facilities.end(); ++facility_it)
        {
            for (int client = 0; client < num_clients; ++client)
            {
                /*  Add this client to the cluster if it is fractionally served by
                the current facility and has not yet been added to any clusters */
                int pos = edge_pos_primal(client, *facility_it, num_clients, num_facilities);
                if (solution_primal[pos] > 0 &&
                    client_assignments.find(client) == client_assignments.end())
                {
                    cluster.clients.insert(client);
                }
            }
        }

        int smallest_cost = INT_MAX, chosen_facility;
        
        /*  Choose a facility to open. In the standard heuristic, the chosen
            facility is the cheapest to open. In the modified heuristic, the
            chosen facility is the one that provides the smallest contribution
            to the objective function */
if (modified)   // Modified heuristic
{
    double smallest_contribution = INT_MAX;

    auto lambda_cmp_contribution = [&](int f1, int f2) -> bool
    {
        return objective_fn_contribution(problem, f1, cluster.clients) <
            objective_fn_contribution(problem, f2, cluster.clients);
    };

    set<int>::const_iterator position = min_element(cluster.facilities.begin(),
        cluster.facilities.end(), lambda_cmp_contribution);

    chosen_facility = *position;
}
else    // Standard heuristic
{
    // Find the facility in the cluster with the smallest opening cost
    for (set<int>::const_iterator facility_it = cluster.facilities.begin();
        facility_it != cluster.facilities.end(); ++facility_it)
    {
        if (problem.facility_costs[*facility_it] < smallest_cost)
        {
            smallest_cost = problem.facility_costs[*facility_it];
            chosen_facility = *facility_it;
        }
    }
}

// Update the mappings
facility_assignments[chosen_facility].push_back(client_seed);
client_assignments[client_seed] = chosen_facility;

// Assign the clients in the cluster to the opened facility
for (set<int>::const_iterator client_it = cluster.clients.begin();
    client_it != cluster.clients.end(); ++client_it)
{
    // Assign the client if there is service from the facility to it
    if (problem.distance_cost[chosen_facility][*client_it] != artificial_inf &&
        *client_it != client_seed)
    {
        client_assignments[*client_it] = chosen_facility;
        facility_assignments[chosen_facility].push_back(*client_it);
    }
}

// Update the cost of the solution
for (set<int>::const_iterator client_it = cluster.clients.begin();
    client_it != cluster.clients.end(); ++client_it)
{
    int distance = problem.distance_cost.at(chosen_facility).at(*client_it);
    if (distance == artificial_inf) distance = 0;
    cost += distance * problem.client_demands[*client_it];
}
cost += problem.facility_costs[chosen_facility];

client_seed = find_cluster_seed(client_assignments, solution_dual, num_clients);
    }
    // ======================================================

    min_cost = cost;
}

// ============================================================================

string print_input_format(const Graph& problem)
{
    int num_facilities = problem.facilities();
    int num_clients = problem.clients();

    stringstream ss;
    ss << num_facilities << " " << num_clients << endl;

    for (int facility = 0; facility < num_facilities; ++facility)
        ss << problem.facility_costs[facility] <<
        ((facility == num_facilities - 1) ? '\n' : ' ');

    for (int client = 0; client < num_clients; ++client)
        ss << problem.client_demands[client] <<
        ((client == num_clients - 1) ? '\n' : ' ');

    for (int facility = 0; facility < num_facilities; ++facility)
        for (int client = 0; client < num_clients; ++client)
            ss << problem.distance_cost[facility][client] <<
            ((client == num_clients - 1) ? '\n' : ' ');

    return ss.str();
}

string print_output_format(int instance, const string& algorithm, 
    int num_facilities, const map<int, int> client_assignments,
    const map<int, vector<int> >& facility_assignments, double cost)
{
    stringstream ss;

    ss << "INSTÂNCIA " << instance << " - " << algorithm << ": Custo " <<
        cost << endl << "Facilities abertas:";

    for (map<int, vector<int> >::const_iterator it = facility_assignments.begin();
        it != facility_assignments.end(); ++it)
    {
        ss << ' ' << (it->first + 1);
    }
    ss << endl;

    for (map<int, vector<int> >::const_iterator it = facility_assignments.begin();
        it != facility_assignments.end(); ++it)
    {
        ss << "Facility f" << (it->first + 1) << " atende clientes:";

        for (vector<int>::size_type i = 0; i < it->second.size(); ++i)
            ss << ' ' << ((it->second)[i] + 1 + num_facilities);

        ss << endl;
    }

    return ss.str();
}

void read_instance(Graph& problem, istream& in, ostream& out)
{
    int num_facilities, num_clients;

    if ( !(in >> num_facilities >> num_clients) )
        throw runtime_error("Could not read from file");

    int buffer;
    for (int facility = 0; facility < num_facilities; ++facility)
    {
        if ( !(in >> buffer) ) throw runtime_error("Could not read from file");
        problem.facility_costs.push_back(buffer);
    }

    for (int client = 0; client < num_clients; ++client)
    {
        if (!(in >> buffer)) throw runtime_error("Could not read from file");
        problem.client_demands.push_back(buffer);
    }

    for (int facility = 0; facility < num_facilities; ++facility)
    {
        problem.distance_cost.push_back( vector<int>(num_clients, 0) );

        for (int client = 0; client < num_clients; ++client)
        {
            if ( !(in >> problem.distance_cost[facility][client]) )
                throw runtime_error("Could not read from file");
        }
    }
}

void read_input(istream& in, ostream& out, bool find_cases = false)
{
    int num_instances;

    if ( !(in >> num_instances) ) 
        throw runtime_error("Could not read number of instances");

    bool standard_found = false, modified_found = false;

    double h_mean = 0, mh_mean = 0, h_worst = INT_MIN, mh_worst = INT_MIN;
    double h_success = 0, mh_success = 0;

    stringstream ss;


    for (int instance = 0; instance < num_instances; ++instance)
    {
        Graph problem;
        read_instance(problem, in, ss);

        map<int, vector<int> > facility_assignment;
        map<int, int> client_assignment;
        double opt_cost, h_cost, mh_cost;

        optimal(problem, client_assignment, facility_assignment, opt_cost);
        ss << print_output_format(instance+1, "Ótima", problem.facilities(), 
            client_assignment, facility_assignment, opt_cost) << endl;

        lp_rounding(problem, client_assignment, facility_assignment, h_cost);
        ss << print_output_format(instance+1, "Heurística", problem.facilities(),
            client_assignment, facility_assignment, h_cost) << endl;

        lp_rounding(problem, client_assignment, facility_assignment, mh_cost, true);
        ss << print_output_format(instance+1, "Heurística Melhorada", 
            problem.facilities(), client_assignment, facility_assignment, mh_cost) << endl;

        ss << "========================================" << endl << endl;

        if (find_cases && h_cost < mh_cost && !standard_found)
        {
            standard_found = true;
            cout << "*** Standard heuristic is better in this instance ***" << endl;
            cout << print_output_format(instance + 1, "Heurística", problem.facilities(),
                client_assignment, facility_assignment, h_cost) << endl;
            cout << print_input_format(problem) << endl;
        }
        else if (find_cases && mh_cost < h_cost && !modified_found)
        {
            modified_found = true;
            cout << "*** Modified heuristic is better in this instance ***" << endl;
            cout << print_output_format(instance + 1, "Heurística Melhorada",
                problem.facilities(), client_assignment, facility_assignment, mh_cost) << endl;
            cout << print_input_format(problem) << endl;
        }

        if (h_cost == opt_cost) h_success++;
        if (mh_cost == opt_cost) mh_success++;

        double h_ratio = h_cost / opt_cost;
        double mh_ratio = mh_cost / opt_cost;

        h_mean += h_ratio;
        mh_mean += mh_ratio;

        h_worst = std::max(h_worst, h_ratio);
        mh_worst = std::max(mh_worst, mh_ratio);
    }

    out << "================================" << endl;
    out << "========== STATISTICS ==========" << endl;
    out << "================================" << endl << endl;
    out << "Heuristic mean ratio: " << (h_mean / num_instances) << endl;
    out << "Modified heuristic mean ratio: " << (mh_mean / num_instances) << endl;
    out << "Heuristic equals OPT: " << (100*(h_success / num_instances)) << "%" << endl;
    out << "Modified heuristic equals OPT: " << (100*(mh_success / num_instances)) << "%" << endl;
    out << "Worst heuristic ratio: " << h_worst << endl;
    out << "Worst modified heuristic ratio: " << mh_worst << endl << endl;

    out << "================================" << endl;
    out << "========== INSTANCES ===========" << endl;
    out << "================================" << endl << endl << ss.str();

}

void read_input(const string& input_filename, const string& out_filename, bool find_cases = false)
{
    ifstream input_file(input_filename);
    ofstream output_file(out_filename);

    if (!input_file) throw runtime_error("Could not open input file");
    if (!output_file) throw runtime_error("Could not open output file");

    read_input(input_file, output_file, find_cases);

    input_file.close();
    output_file.close();
}

void random_problems(int problems_per_configuration, ostream& out)
{

    out << "=====================================" << endl;
    out << "========== Configuration 1 ==========" << endl;
    out << "=====================================" << endl << endl;
    for (int i = 0; i < problems_per_configuration; i++)
        out << print_input_format( random_graph(4, 10) ) << endl;

    out << "=====================================" << endl;
    out << "========== Configuration 2 ==========" << endl;
    out << "=====================================" << endl << endl;
    for (int i = 0; i < problems_per_configuration; i++)
        out << print_input_format( random_graph(5, 12) ) << endl;

    out << "=====================================" << endl;
    out << "========== Configuration 3 ==========" << endl;
    out << "=====================================" << endl << endl;
    for (int i = 0; i < problems_per_configuration; i++)
        out << print_input_format( random_graph(6, 14) ) << endl;

    out << "=====================================" << endl;
    out << "========== Configuration 4 ==========" << endl;
    out << "=====================================" << endl << endl;
    for (int i = 0; i < problems_per_configuration; i++)
        out << print_input_format( random_graph(7, 16) ) << endl;
}

int main()
{
    env.set(GRB_IntParam_OutputFlag, 0);
    
    try
    {
        string choice;
        cout << "Type \"f\" for file inputs, or \"r\" for random problems:" << endl;
        
        do { std::cin >> choice; } while (choice != "r" && choice != "f");

        if (choice == "f")
            read_input("Projeto_UFL_input.txt", "Projeto_UFL_output_jgsma.txt", true);
        else
        {
            ofstream output_file("Projeto_UFL_jgsma_400_inputs.txt");
            random_problems(400, output_file);
            output_file.close();
        }
    }
    catch (std::exception& e)
    {
        cout << "ERROR => " << e.what() << endl;
    }
}

#endif
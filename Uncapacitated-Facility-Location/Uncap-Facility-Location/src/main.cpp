#ifndef __GUARD_MAIN__
#define __GUARD_MAIN__

#include <algorithm>
#include <climits>
#include <iterator>
#include <map>
#include <numeric>
#include <vector>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>

#include "pilal.h"
#include "simplex.h"

using std::accumulate;
using std::back_inserter;
using std::copy;
using std::map;
using std::min_element;
using std::ostream;
using std::pair;
using std::runtime_error;
using std::set;
using std::string;
using std::stringstream;
using std::transform;
using std::vector;

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

    void print_constraints(ostream& strm, const char constraint_type,
        const vector<vector<double> >& constraints) const;
    
    void to_lib_syntax(ostream& strm, bool do_maximize, bool relaxed,
        const vector<double>& objective_function,
        const vector<vector<double> >& equality_constraints,
        const vector<vector<double> >& lesser_than_constraints,
        const vector<vector<double> >& greater_than_constraints) const;
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
    solution.clear();

    vector<vector<double> > eq_constraints, lt_constraints, gt_constraints;

    for (vector<vector<double> >::size_type i = 0; i < constraints.size(); ++i)
    {
        switch (constraint_types[i])
        {
        case '=': eq_constraints.push_back(constraints[i]); break;
        case '<': lt_constraints.push_back(constraints[i]); break;
        case '>': gt_constraints.push_back(constraints[i]); break;
        }
    }

    stringstream ss;

    to_lib_syntax(ss, do_maximize, true, objective_vector, eq_constraints, 
        lt_constraints, gt_constraints);

    optimization::Simplex solver("LP Problem");
    solver.load_problem(ss);

    solver.solve();

    if ( solver.must_be_fixed() )
        throw runtime_error("Problem formulation is incorrect");

    if ( solver.has_solutions() )
    {
        if ( !solver.is_unlimited() ) 
            // Return the solution tht was found
            solver.get_solution(solution, value);
        else 
            throw runtime_error("Problem is unlimited");
    }
    else
        throw runtime_error("Problem is overconstrained");
}

void Linear_Program::print_constraints(ostream& strm, const char constraint_type,
    const vector<vector<double> >& constraints) const
{
    for (int row = 0; row < constraints.size(); ++row)
    {
        int size = (int) constraints[row].size() - 1;
        for (int col = 0; col < size; ++col)
            strm << constraints[row][col] << '\t';

        strm << " " << constraint_type << " " << constraints[row].back()
            << std::endl;
    }
}

void Linear_Program::to_lib_syntax(ostream& strm, bool do_maximize, bool relaxed,
    const vector<double>& objective_function,
    const vector<vector<double> >& equality_constraints,
    const vector<vector<double> >& lesser_than_constraints,
    const vector<vector<double> >& greater_than_constraints) const
{
    string limit = relaxed ? "inf" : "1";
    string min_max = do_maximize ? "maximize" : "minimize";
    
    strm << "[METADATA]" << std::endl;
    strm << "vars " << objective_function.size() << std::endl;
    strm << "[VARIABLES]" << std::endl;

    for (int var_num = 1; var_num <= objective_function.size(); ++var_num)
        strm << "0    var" << var_num << "    " << limit << std::endl;

    strm << "[CONSTRAINTS]" << std::endl;

    print_constraints(strm, '=', equality_constraints);
    print_constraints(strm, '<', lesser_than_constraints);
    print_constraints(strm, '>', greater_than_constraints);

    strm << "[OBJECTIVE]" << std::endl;
    strm << min_max;

    accumulate(objective_function.begin(), objective_function.end(), &strm,
        [](ostream* stream, double d) -> ostream* { return &(*stream << " " << d); });
    strm << std::endl;
}
// ============================================

// ============================================================================



// ============================================================================
// ================================== GRAPH ===================================
// ============================================================================

struct Graph {
    map<int, map<int, int> > distance_cost;
    vector<int> facility_costs;
    vector<int> client_demands;

    int facilities() const { return (int)facility_costs.size(); }

    int clients() const { return (int)client_demands.size(); }

private:
    static int get_key(const pair<int, int>& map_it) { return map_it.first; }
};

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

void optimal(const Graph& problem, vector<double>& open_facilities,
    map<int, int>& client_assignments, double& min_cost)
{
    open_facilities.clear();
    client_assignments.clear();

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
    
    for (int i = 0; i < num_facilities; ++i)
        open_facilities.push_back( ((1ll << i) & best_opened) > 0 );
    
    for (int client = 0; client < num_clients; ++client)
    {
        for (int facility = 0; facility < num_facilities; ++facility)
        {
            int position = edge_pos_primal(client, facility, num_clients);
            if ( best_solution[position] == 1 )
            {
                client_assignments[client] = facility;
                break;
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
        return acc + problem.distance_cost.at(facility).at(client) * 
            problem.client_demands[client]; 
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
            constraint[client] = 1.;

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

void lp_rounding(const Graph& problem, vector<double>& open_facilities,
    map<int, int>& client_assignments, double& min_cost, bool modified = false)
{
    open_facilities.clear();
    client_assignments.clear();

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

    open_facilities.resize(problem.facilities());

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

        // Mark the cheapest facility as opened
        open_facilities[chosen_facility] = 1;

        // Assign all clients in the cluster to the opened facility
        for (set<int>::const_iterator client_it = cluster.clients.begin();
            client_it != cluster.clients.end(); ++client_it)
        {
            client_assignments[*client_it] = chosen_facility;
        }

        // Update the cost of the solution
        for (set<int>::const_iterator client_it = cluster.clients.begin();
            client_it != cluster.clients.end(); ++client_it)
        {
            cost += problem.distance_cost.at(chosen_facility).at(*client_it) *
                problem.client_demands[*client_it];
        }
        cost += problem.facility_costs[chosen_facility];

        client_seed = find_cluster_seed(client_assignments, solution_dual, num_clients);
    }
    // ======================================================

    min_cost = cost;
}

// ============================================================================



void print_solution(const vector<double>& open_facilities,
    const map<int, int>& client_assignment, double cost)
{
    std::cout << "Total Cost: " << cost << std::endl << std::endl;

    std::cout << "Solution: " << std::endl;
    for (int i = 0; i < open_facilities.size(); ++i)
        std::cout << "\tFacility #" << (i + 1) << "\t" <<
        open_facilities[i] << std::endl;
    std::cout << std::endl;

    std::cout << "Client assignment: " << std::endl;
    for (map<int, int>::const_iterator it = client_assignment.begin();
        it != client_assignment.end(); ++it)
    {
        std::cout << "\tClient #" << (it->first + 1) << "\t" <<
            "Facility #" << (it->second + 1) << std::endl;
    }
}

int main()
{
    Graph g;
    g.distance_cost[0][0] = 1;
    g.distance_cost[0][1] = 10;
    g.distance_cost[0][2] = 6;
    g.distance_cost[0][3] = 6;
    g.distance_cost[1][0] = 7;
    g.distance_cost[1][1] = 4;
    g.distance_cost[1][2] = 2;
    g.distance_cost[1][3] = 9;
    g.distance_cost[2][0] = 7;
    g.distance_cost[2][1] = 9;
    g.distance_cost[2][2] = 10;
    g.distance_cost[2][3] = 10;
    g.distance_cost[3][0] = 2;
    g.distance_cost[3][1] = 6;
    g.distance_cost[3][2] = 3;
    g.distance_cost[3][3] = 10;
    g.distance_cost[4][0] = 7;
    g.distance_cost[4][1] = 9;
    g.distance_cost[4][2] = 6;
    g.distance_cost[4][3] = 10;

    g.facility_costs = { 7, 10, 6, 3, 10 };
    g.client_demands = { 10, 3, 4, 6 };

    vector<double> open_facilities;
    map<int, int> client_assignment;
    double cost;

    optimal(g, open_facilities, client_assignment, cost);

    std::cout << "===== OPTIMAL =====" << std::endl;
    print_solution(open_facilities, client_assignment, cost);
    std::cout << std::endl << std::endl;

    lp_rounding(g, open_facilities, client_assignment, cost);
    std::cout << "===== LP ROUNDING =====" << std::endl;
    print_solution(open_facilities, client_assignment, cost);
    std::cout << std::endl << std::endl;

    lp_rounding(g, open_facilities, client_assignment, cost, true);
    std::cout << "===== MODIFIED LP ROUNDING =====" << std::endl;
    print_solution(open_facilities, client_assignment, cost);
}

#endif
#ifndef __GUARD_MAIN__
#define __GUARD_MAIN__

#include <algorithm>
#include <climits>
#include <iterator>
#include <map>
#include <numeric>
#include <vector>
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
using std::ostream;
using std::pair;
using std::runtime_error;
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
// =============================== BRUTE-FORCE ================================
// ============================================================================

inline int edge_pos(int client, int facility, int num_clients)
{
    return (facility * num_clients) + client;
}

inline bool is_open(unsigned long long configuration, int facility)
{
    return ((1ll << facility) & configuration) > 0;
}

void optimal(const Graph& problem, vector<bool>& open_facilities,
    map<int, int>& client_assigned_facility, double& min_cost)
{
    open_facilities.clear();
    client_assigned_facility.clear();

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
                int pos = edge_pos(client, facility, num_clients);
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
                int position = edge_pos(client, facility, num_clients);
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

                int pos = edge_pos(client, facility, num_clients);
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
            int position = edge_pos(client, facility, num_clients);
            if ( best_solution[position] > 0.5 )
            {
                client_assigned_facility[client] = facility;
                break;
            }
        }
    }

}


// ============================================================================


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

    vector<bool> open_facilities;
    map<int, int> client_assignment;
    double max;

    optimal(g, open_facilities, client_assignment, max);

}

#endif
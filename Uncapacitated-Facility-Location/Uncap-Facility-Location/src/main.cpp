#ifndef __GUARD_MAIN__
#define __GUARD_MAIN__

#include <vector>
#include <stdexcept>
#include <string>

using std::vector;
using std::runtime_error;
using std::string;


// ============================================================================
// ============================== LINEAR PROGRAM ==============================
// ============================================================================

// ============= CLASS DECLARATION =============
class Linear_Program {
public:    
    Linear_Program(const vector<double>& objective_vector, const vector<vector<double> >& a_matrix,
        const vector<double>& b_vector) : objective_vector(objective_vector), a_matrix(a_matrix),
        b_vector(b_vector) { validate_values(&objective_vector, &a_matrix, &b_vector); }

    Linear_Program(const vector<double>& objective, const vector<vector<double> >& constraints);

    void solve(vector<double>& solution, double& max) const;

private:
    vector<double> objective_vector;    // Coefficients of the objective function to be maximized
    vector<vector<double> > a_matrix;   // A matrix of the linear program in standard form
    vector<double> b_vector;            // b vector of the linear program in standard form

    void validate_values(const vector<double>* objective_vector,
        const vector<vector<double> >* a_matrix, const vector<double>* b_vector = nullptr) const;

    void set_simplex(vector<vector<double>>& simplex) const;
    
    bool get_pivots(vector<vector<double>> simplex, int& pivot_col, 
        int& pivot_row, bool& no_solution) const;

    void do_simplex(vector<vector<double>> simplex, vector<double>& x, double& max) const;
};
// =============================================


// ============= CLASS DEFINITION =============
Linear_Program::Linear_Program(const vector<double>& objective_vector,
    const vector<vector<double> >& constraints): objective_vector(objective_vector)
{
    validate_values(&objective_vector, &constraints);

    for (vector<vector<double> >::const_iterator it = constraints.begin();
        it != constraints.end(); ++it)
    {
        const vector<double>& constraint = *it;
        
        // Add the b value contained in the end of this constraint
        double b_constant = constraint.back();
        b_vector.push_back(b_constant);
        
        // Add the coefficients in this constraint
        a_matrix.push_back( 
            vector<double>(constraint.begin(), --constraint.end()) );
    }
}

void Linear_Program::validate_values(const vector<double>* objective_vector,
    const vector<vector<double> >* a_matrix, 
    const vector<double>*  b_vector) const
{
    const string& error_message = "Invalid Linear Programming configuration";

    // Check if either vector is empty
    if (objective_vector->empty() || a_matrix->empty() ||
        (b_vector != nullptr && b_vector->empty()))
    {
        throw runtime_error(error_message);
    }

    // Check if objective vector of coefficients contains coefficients 
    // for all variables
    int offset = (b_vector == nullptr) ? 1 : 0;
    if (objective_vector->size() != a_matrix->at(0).size() - offset)
        throw runtime_error(error_message);

    // Check if there is one b value for each equation
    if (b_vector != nullptr && b_vector->size() != a_matrix->size())
        throw runtime_error(error_message);

    // Check if all constraints contain coefficients for all variables
    vector<double>::size_type num_variables = a_matrix->at(0).size();
    for (vector<double>::size_type i = 1; i < a_matrix->size(); ++i)
    {
        if (a_matrix->at(i).size() != num_variables)
            throw runtime_error(error_message);
    }
}

void Linear_Program::set_simplex(vector<vector<double>>& simplex) const
{
    simplex.clear();

    int num_variables = (int) objective_vector.size();
    int num_equations = (int) a_matrix.size();
    int num_cols = num_variables + num_equations + 1 + 1;


    for (int irow = 0; irow < num_equations; irow++)
    {
        vector<double> row(num_cols, 0);
        for (int icol = 0; icol < num_variables; icol++)
        {
            row[icol] = a_matrix[irow][icol];
        }
        row[num_variables + irow] = 1;
        row[num_cols - 1] = b_vector[irow];


        simplex.push_back(row);
    }


    vector<double> last_row(num_cols, 0);
    for (int i_var = 0; i_var < num_variables; i_var++)
    {
        last_row[i_var] = 0 - objective_vector[i_var];
    }
    last_row[num_variables + num_equations] = 1;
    simplex.push_back(last_row);
}

bool Linear_Program::get_pivots(vector<vector<double>> simplex, int& pivot_col,
    int& pivot_row, bool& no_solution) const
{
    int num_rows = (int) simplex.size();
    int num_cols = (int) simplex[0].size();
    int num_variables = num_cols - num_rows - 1;

    no_solution = false;

    double min = 0;
    for (int icol = 0; icol < num_cols - 2; icol++)
    {
        double value = simplex[num_rows - 1][icol];
        if (value < min)
        {
            pivot_col = icol;
            min = value;
        }
    }


    if (min == 0)
        return false;


    double min_ratio = 0.0;
    bool first = true;
    for (int irow = 0; irow < num_rows - 1; irow++)
    {
        double value = simplex[irow][pivot_col];

        if (value > 0.0)
        {
            double col_value = simplex[irow][num_cols - 1];
            double ratio = col_value / value;


            if ((first || ratio < min_ratio) && ratio >= 0.0)
            {
                min_ratio = ratio;
                pivot_row = irow;
                first = false;
            }
        }
    }


    no_solution = first;
    return !no_solution;
}

void Linear_Program::do_simplex(vector<vector<double>> simplex,
    vector<double>& x, double& max) const
{
    x.clear();

    int pivot_col, pivot_row;
    int num_rows = (int) simplex.size();
    int num_cols = (int) simplex[0].size();

    bool no_solution = false;
    while (get_pivots(simplex, pivot_col, pivot_row, no_solution))
    {
        double pivot = simplex[pivot_row][pivot_col];

        for (int icol = 0; icol < num_cols; icol++)
        {
            simplex[pivot_row][icol] /= pivot;
        }

        for (int irow = 0; irow < num_rows; irow++)
        {
            if (irow != pivot_row)
            {
                double ratio = -1 * simplex[irow][pivot_col];
                for (int icol = 0; icol < num_cols; icol++)
                {
                    simplex[irow][icol] = simplex[pivot_row][icol] * ratio + simplex[irow][icol];
                }
            }
        }
    }

    if (no_solution)
        return;
    
    //optimo!!!
    max = simplex[num_rows - 1][num_cols - 1];
    int num_variables = num_cols - num_rows - 1;
    x = vector<double>(num_variables, 0);

    for (int icol = 0; icol < num_variables; icol++)
    {
        bool is_unit = true;
        bool first = true;
        double value;
        for (int j = 0; j < num_rows; j++)
        {
            if (simplex[j][icol] == 1.0 && first)
            {
                first = false;
                value = simplex[j][num_cols - 1];
            }
            else if (simplex[j][icol] != 0.0)
            {
                is_unit = false;
            }
        }

        if (is_unit && !first)
            x[icol] = value;
        else
            x[icol] = 0.0;
    }
}

void Linear_Program::solve(vector<double>& solution, double& max) const
{
    solution.clear();

    vector<vector<double> > simplex;

    set_simplex(simplex);
    do_simplex(simplex, solution, max);
}
// ============================================

// ============================================================================



// ============================================================================
// ================================== GRAPH ===================================
// ============================================================================



// ============================================================================


// Test case
void run_test()
{
    vector<double> maxFunc = { 0.5 , 3.0 , 1.0 , 4.0 };

    vector<vector<double>> A = 
    {
        {1.0, 1.0, 1.0, 1.0, 40}, 
        {-2.0, -1.0, 1.0, 1.0, 10}, 
        {0.0, 1.0, 0.0, -1.0, 10} 
    };
    
    Linear_Program lp(maxFunc, A);

    double max;
    vector<double> result;

    lp.solve(result, max);

    printf("Result: Max = %f\n", max);
    for (int i = 0; i < result.size(); i++)
    {
        printf("x%d = %f ; ", i, result[i]);
    }
    printf("\n----------------------\n");
}

int main() { run_test(); }

#endif
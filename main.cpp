#include <fstream>
#include <iostream>
#include <string>
#include <vector>

using std::cout;
using std::endl;
using std::fstream;
using std::getline;
using std::istream;
using std::string;
using std::vector;

struct Edge;

struct Vertex {
    unsigned int id;   
    vector<Edge*> edges; 
};

struct Edge {
    Vertex* orig;
    Vertex* dest;
    unsigned int capacity;
};

istream& read_problem_instance(istream& stream, vector<Vertex*>& vertices, 
      vector<Edge*>& edges) {
    int num_vertices, num_edges;

    if (!stream)
        return stream;

    if (!(stream >> num_edges) || !(stream >> num_vertices))
        cout << "ERROR READING INPUT DATA";
    
    bool marked_vertices[num_vertices];

    for (int i = 0; i < num_edges; i++) {
        unsigned int origin, destination, capacity;
        
        if (!(stream >> origin) || !(stream >> destination) || 
              !(stream >> capacity))
            cout << "ERROR READING INPUT DATA";

        Vertex* orig_vertex = new Vertex;
        Vertex* dest_vertex = new Vertex;
        Edge* edge = new Edge;
        
        dest_vertex->id = destination - 1;
        orig_vertex->id = origin - 1;
        orig_vertex->edges.push_back(edge);
        edge->orig = orig_vertex;
        edge->dest = dest_vertex;
        edge->capacity = capacity;

        edges.push_back(edge);
        // Add both origin and destination vertices if they are not repetitions
        if ( !marked_vertices[orig_vertex->id] )
            vertices.push_back(orig_vertex);
        if ( !marked_vertices[dest_vertex->id] )
            vertices.push_back(dest_vertex);

        // Set both vertices as marked
        marked_vertices[orig_vertex->id] = true;
        marked_vertices[dest_vertex->id] = true;
    }

    return stream;
}

int main() {
    fstream input_file;
    input_file.open("in.txt", fstream::in);

    if (!input_file)
        cout << "COULD NOT OPEN FILE" << endl;

    vector<Vertex*> vertices;
    vector<Edge*> edges;

    read_problem_instance(input_file, vertices, edges);

    return 0;
}
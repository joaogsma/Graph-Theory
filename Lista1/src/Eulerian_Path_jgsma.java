import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.LinkedList;

public class Eulerian_Path_jgsma {
	
    // ========================================================================
	// ============================ SINGLETON CODE ============================
	// ========================================================================
    private static Eulerian_Path_jgsma instance; 
	
	private Eulerian_Path_jgsma(){}
	
	
	public static Eulerian_Path_jgsma get_instance(){
		if (instance == null)
			instance = new Eulerian_Path_jgsma();
		
		return instance;
	}
	
	
	private Graph create_graph_instance() {
		return new Graph();
	}
	
	
	public static Graph create_graph() {
		return get_instance().create_graph_instance();
	}
	// ========================================================================
	
	
	
    // ========================================================================
	// ========================== UTILITY FUNCTIONS ===========================
	// ========================================================================
    // Converts a vertex list to string format, e.g.: 1 -> 2 -> 3 -> 4
    public static String vertex_list_to_str(List<Vertex> list) {
    	StringBuilder sb = new StringBuilder();
    	int pos = 0;
    	
    	for (Vertex v : list) {
    		sb.append(v.id);
    		if (pos != list.size() - 1) { sb.append(" -> "); }
    		pos++;
    	}
    	
    	return sb.toString();
    }
	
    
    /* Reads a text input file containing a graph. In the file, I am
     * assuming that in the line corresponding to a given vertex, its
     * adjacent vertices are separated by spaces. I am also assuming
     * that vertex numbering starts at 1. */
    public static Graph load_graph(String filename) throws IOException {
    	Graph g = create_graph();
    	
    	File f = new File(filename);
    	FileReader fr = new FileReader(f);
    	BufferedReader br = new BufferedReader(fr);
    	
    	// Read the number of vertices and create them
    	int n_vertices = Integer.parseInt(br.readLine());
    	for (int id = 1; id <= n_vertices; id++)
    		g.add_vertex(id);
    	
    	// Read the edges and add them
    	for (int id = 1; id <= n_vertices; id++) {
    		// Read the adjacent vertices
    		String line = br.readLine();
    		String[] adjacents = line.split(" ");
    		
    		// Create the edges
    		for (String str : adjacents) {
    			int id_ = Integer.parseInt(str);
    			
    			// Check if this is a new edge
    			boolean new_edge = true;
    			for (Edge e : g.vertices.get(id).edges) {
    				if (e.v1.id == id_ || e.v2.id == id_) {
    					new_edge = false;
    					break;
    				}
    			}
    			
    			// Add the edge if its a new one
    			if (new_edge)
    				g.add_edge(id, id_);
    		}
    	}
    	
    	br.close();
    	return g;
    }
    // ========================================================================
	
    
	
    // ========================================================================
	// ============================ GRAPH CLASSES =============================
	// ========================================================================
    class Vertex {
        List<Edge> edges;
        final int id;

        public Vertex(int id) {
            this.id = id;
            edges = new LinkedList<Edge>();
        }

        
        // Returns the vertex adjacent to this, with respect to a given edge
        public Vertex adjacent(Edge e) {
            if (e.v1 == this)
                return e.v2;

            if (e.v2 == this)
                return e.v1;

            throw new RuntimeException("Edge not incident on specified vertex");
        }
    }

	
    class Edge {
        Vertex v1, v2;

        public Edge(Vertex v1, Vertex v2) {
            this.v1 = v1;
            this.v2 = v2;
        }
    }

    
    class Graph {
        HashMap<Integer, Vertex> vertices;
        HashSet<Edge> marked;
        
        /* This class represents a simple linked list, but with each 
         * node having access to the end of the list. The reason for 
         * this is to allow constant time addition to the end of the
         * list, as is done frequently in the computation of the 
         * euler circuits, when concatenating linked lists. */
        class List_Node<E> {
            E data;
            List_Node<E> next, tail;

            List_Node(E v) {
                data = v;
                tail = this;
            }
            
            void append(List_Node<E> node) {
            	tail.next = node;
            	tail = node.tail;
            }
        }

        public Graph() {
            vertices = new HashMap<Integer, Vertex>();
            marked = new HashSet<Edge>();
        }

        
        // Adds a vertex to the graph
        public void add_vertex(int id) {
        	vertices.put(id, new Vertex(id));
        }

        
        public void add_edge(int id1, int id2) {
        	// Check for invalid vertices
        	if (!vertices.containsKey(id1) || !vertices.containsKey(id2))
        		throw new RuntimeException("Graph does not contain vertex");
        	
        	// Get corresponding vertices
        	Vertex v1 = vertices.get(id1);
        	Vertex v2 = vertices.get(id2);
        	
        	// Create edge and link vertices to it
        	Edge e = new Edge(v1, v2);
        	v1.edges.add(e);
        	v2.edges.add(e);
        }

        
        // Computes an euler circuit starting at the vertex 1
        public List<Vertex> euler_circuit_hierholzer() {
        	return euler_circuit_hierholzer(1);
        }
        
        
        // Computes an euler circuit starting at a given vertex
        public List<Vertex> euler_circuit_hierholzer(int start_vertex_id) {
        	if (vertices.isEmpty())
        		throw new RuntimeException("Empty graph");
        	
        	marked.clear(); 	// Clean-up phase
        	List_Node<Vertex> circuit_head = rec_euler_circuit_hierholzer(
        			vertices.get(start_vertex_id));
        	
            List<Vertex> complete_circuit = new LinkedList<Vertex>();
            while (circuit_head != null) {
                complete_circuit.add(circuit_head.data);
                circuit_head = circuit_head.next;
            }

        	return complete_circuit;
        }
        
        
        /* This function is the actual computation of Hierholzen's algorithm. It
         * recursively computes a cycle C in a given connected graph G (this 
         * Graph instance with a subset of its edges) and merges it with the 
         * euler circuits of the connected components of G - {edges of C}. */
        private List_Node<Vertex> rec_euler_circuit_hierholzer(Vertex v) {
        	List_Node<Vertex> euler_circuit = null;

        	// Check if vertex has odd degree
    		if (v.edges.size() % 2 != 0)
    			throw new RuntimeException("Vertex with odd degree");
        	
    		// Find a cycle starting at the vertex v
        	List_Node<Vertex> cycle = find_cycle(v);
        	        	
        	// Stop condition. Return only the vertex if there are no cycles
        	if (cycle.next == null)
        		return cycle;
        	
        	// Join cycles of connected components
        	List_Node<Vertex> head = cycle;
        	while (head != null) {
        		// Compute the subgraph's euler circuit
        		List_Node<Vertex> circuit = rec_euler_circuit_hierholzer(head.data);
                
        		/* If the circuit has more than one vertex, close it by adding the 
        		 * initial vertex at the end */
        		if (circuit.next != null) 
        			circuit.append(new List_Node<Vertex>(head.data)); 
        		
        		// Add the subgraph's euler circuit to the complete circuit
                if (euler_circuit == null)
                	euler_circuit = circuit;
                else
                	euler_circuit.append(circuit);
                
                head = head.next;
        	}
        	
        	return euler_circuit;
        }
        
        
        // Finds a cycle in the graph
        private List_Node<Vertex> find_cycle(Vertex v) {
        	// Check if the graph contains the vertex 
        	if (!vertices.containsKey(v.id))
        		throw new RuntimeException("Graph does not contain vertex");
        	
        	List_Node<Vertex> cycle = null;
        	
        	// Go through vertices until arrival back at v
        	Vertex current = v;
        	do {
        		// Add current vertex to the cycle
        		if (cycle == null)
        			cycle = new List_Node<Vertex>(current);
        		else 
	        		cycle.append(new List_Node<Vertex>(current));
        		        		
        		// Find an unmarked edge
        		for (Edge e : current.edges) {
        			if (!marked.contains(e)) {
        				// Mark the edge
        				marked.add(e);
        				// Move current to the adjacent vertex
        				current = current.adjacent(e);
        				break;
        			}
        		}
        		
        	} while(current != v);
        	
        	return cycle;
        }
    }
    // ========================================================================
    
    
    
    // ========================================================================
    // ================================= MAIN =================================
    // ========================================================================
    public static void main(String[] args) throws IOException {
    	String filename;
    	Graph g;
    	
    	try {
    		filename = args[args.length - 1];
    		g = load_graph(filename);
    	} catch (ArrayIndexOutOfBoundsException e) {
    		System.out.println("Please specify the input file");
    		return;
    	} catch (IOException e) {
    		if (e instanceof FileNotFoundException)
    			System.out.println("Could not find file");
    		else 
    			System.out.println("Could not read file");
    		
    		return;
    	}
    	
    	List<Vertex> euler_circuit = g.euler_circuit_hierholzer();
    	String euler_circuit_str = vertex_list_to_str(euler_circuit);
    	System.out.println(euler_circuit_str + (" -> " + euler_circuit.get(0).id));
    }
    // ========================================================================

}


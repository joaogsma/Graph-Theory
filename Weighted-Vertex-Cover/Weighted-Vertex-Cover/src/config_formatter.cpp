#include <iostream>
#include <sstream>
#include <string>

using std::string;
using std::cin;
using std::cout;
using std::endl;


/*	The graph drawer used is the one availeable at:
	http://g.ivank.net/ */

int to_graph_drawer_format()
{
	string in;
	string out;

	cout << "Type the graph string followed by a EOF (Ctrl + Z in Windows):" << endl;

	int line = 0;
	while( std::getline(cin, in) )
	{
		line++;
		
		{
			std::stringstream ss;
			ss << (line == 1? "" : ",") << line << "-";
			out += ss.str();
		}

		for (int i = 0; i < in.size(); i++)
		{
			if ( in[i] == ' ' )
			{
				std::stringstream ss;
				ss << "," << line << "-";
				out += ss.str();
			}
			else
			{
				out += in[i];
			}
		}
	}
    
	std::stringstream ss;
	ss << line << ":" << out;

	out = ss.str();

    cout << out << endl;

}
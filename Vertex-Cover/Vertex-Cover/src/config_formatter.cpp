#include <iostream>
#include <sstream>
#include <string>

using std::string;
using std::cin;
using std::cout;
using std::endl;


int main()
{
	string in;
	string out;

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
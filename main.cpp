#include "xtalcomp.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <functional>
#include <cctype>
#include <locale>

bool verbose = false;

// trim from start
static inline std::string &ltrim(std::string &s)
{
    s.erase(s.begin(), std::find_if(s.begin(), s.end(),
            std::not1(std::ptr_fun<int, int>(std::isspace))));
    return s;
}

// trim from end
static inline std::string &rtrim(std::string &s)
{
    s.erase(std::find_if(s.rbegin(), s.rend(),
            std::not1(std::ptr_fun<int, int>(std::isspace))).base(), s.end());
    return s;
}

// trim from both ends
static inline std::string &trim(std::string &s) 
{
    return ltrim(rtrim(s));
}

// Represents a crystal structure
struct structure
{
	XcMatrix cell;
	std::vector<unsigned int> types;
	std::vector<XcVector> positions;
};

// Load a castep .cell file into a structure
structure load_cell(std::string filename)
{
	// Read in file line by line
	std::ifstream file(filename);
	std::string line;
	std::vector<std::string> lines;
	while(std::getline(file, line))
	{
		trim(line);
		if (line.find("#") == 0) continue; // Skip comments
		if (line.find(";") == 0) continue; // Skip comments
		if (line.find("!") == 0) continue; // Skip comments
		if (line.find("COMMENT") == 0) continue; // Skip comments
		std::transform(line.begin(), line.end(), line.begin(), ::tolower); // Make lowercase
		lines.push_back(line);
	}

	// Parse the structure from the cell file
	structure s;

	std::vector<std::string> names;
	for (int i=0; i<lines.size(); ++i)
	{
		std::string line = lines[i];

		// Parse a block
		if (line.find("%block") != std::string::npos)
		{
			// Parse a lattice_cart block into an XcMatrix
			if (line.find("lattice_cart") != std::string::npos)
			{
				
				// Locate the start and end of the lattice
				// by finding the %endblock and working backwards
				// (as the first line of the block can optionally
				//  specify the units)
				int start = -1;
				int end   = -1;
				for (int j=i+1; j<lines.size(); ++j)
					if (lines[j].find("%endblock") != std::string::npos)
					{
						start = j-3;
						end   = j;
						break;
					}

				if (start < 0)
				{
					std::cout << "Error parsing lattice_cart, could not find %endblock statement!";
					throw;
				}

				double lattice[3][3];
				for (int j=start; j<end; ++j)
				{
					std::stringstream stream(lines[j]);
					double x, y, z;
					stream >> x;
					stream >> y;
					stream >> z;

					lattice[j-start][0] = x;
					lattice[j-start][1] = y;
					lattice[j-start][2] = z;
				}

				s.cell = XcMatrix(lattice);
			}

			// Parse a positions_frac block into an XcVector
			else if (line.find("positions_frac") != std::string::npos)
			{
				int j = i+1;
				std::string subline = lines[j];

				// Keep parsing until we hit %endblock
				while(subline.find("%endblock") == std::string::npos)
				{
					std::stringstream stream(subline);
					std::string name;
					double x, y, z;
					stream >> name >> x >> y >> z;

					int type = -1;
					for (int k = 0; k<names.size(); ++k)
						if (names[k].compare(name) == 0)
						{
							type = k;
							break;
						}

					if (type < 0)
					{
						names.push_back(name);
						type = names.size() - 1;
					}

					s.types.push_back(type);
					s.positions.push_back(XcVector(x,y,z));

					// Move to next atom
					++j;
					subline = lines[j];
				}
			}	

			else if (line.find("positions_cart") != std::string::npos)
			{
				std::cout << "Error, positions_cart not supported: email mjh261@cam.ac.uk\n";
				throw;
			}
		}
	}

	if (verbose)
	{
		std::cout << "Cell " << filename << "\n";
		std::cout << "Lattice:\n";
		for (int i=0; i<3; ++i)
			std::cout << s.cell.row(i).x() << " " 
				  << s.cell.row(i).y() << " " 
				  << s.cell.row(i).z() << "\n";

		std::cout << "Atoms (" << s.types.size() << ")\n";
		for (int i=0; i<s.types.size(); ++i)
			std::cout << names[s.types[i]] << " " 
				  << s.positions[i].x() << " "
				  << s.positions[i].y() << " "
				  << s.positions[i].z() << "\n";
		std::cout << "\n";
	}

	return s;
}

int main(int argc, char** argv)
{
	if (argc < 3)
	{
		std::cout << "Error, I require at least two arguments: the cells to compare!\n";
		return -1;
	}

	double cart_tol  = 0.05;
	double angle_tol = 0.25;

	for (int i=3; i<argc; ++i)
	{
		std::string arg = argv[i];
		int ep = arg.find("=");

		if (arg.find("-cart_tol=") != std::string::npos)
			cart_tol = std::stod(arg.substr(ep+1, arg.size()-ep-1));
		
		else if (arg.find("-angle_tol=") != std::string::npos)
			angle_tol = std::stod(arg.substr(ep+1, arg.size()-ep-1));

		else if (arg.find("-v") != std::string::npos)
			verbose = true;
	}

	structure s1 = load_cell(argv[1]);
	structure s2 = load_cell(argv[2]);
	bool same = XtalComp::compare(s1.cell, s1.types, s1.positions,
				      s2.cell, s2.types, s2.positions,
				      NULL, cart_tol, angle_tol, true);

	std::string same_str = "Not equivalent";
	if (same) same_str = "Equivalent";
	std::cout << same_str << " to within cart_tol = " << cart_tol << " angle_tol = " << angle_tol << "\n";

}

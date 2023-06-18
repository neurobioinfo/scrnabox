//============================================================================
// Name        : QueueInterpreter.cpp
// Author      : Alex D L
// Version     :
// Copyright   : copyright notice
// Description : Hello World in C, Ansi-style
//============================================================================


#include <stdio.h>
#include <stdlib.h>
#include <fstream>
#include <ctime>

#include <boost/program_options.hpp>
using namespace boost::program_options;

#include <boost/foreach.hpp>
#include <boost/tokenizer.hpp>


#include <iostream>
#include <iterator>
using namespace std;

// A helper function to simplify the main part.
template<class T>
ostream& operator<<(ostream& os, const vector<T>& v)
{
	copy(v.begin(), v.end(), ostream_iterator<T>(cout, " "));
	return os;
}

std::string time_t_to_string(time_t t)
{
    std::stringstream sstr;
    sstr << t;
    return sstr.str();
}

int main(int ac, char* av[]) {

	int returnVal=EXIT_FAILURE;

	try {
		string input_file, account, name;
		vector <string> resources, variables, dependencies;
		bool exportAll, verbose;

		string time_milli = time_t_to_string(time(0));

		// Declare a group of options that will be allowed only on command line
		options_description generic("Generic options");
		generic.add_options()
            		("version", "print version string")
            		("help", "produce help message")
            		("verbose", value<bool>(&verbose)->implicit_value(true)->default_value(false), "Show verbose output")
            		;

		// Declare a group of options that will be allowed both on command line and in config file
		options_description config("Configuration");
		config.add_options()
		  			  ("account,A", value<string>(&account), "Account name")
		  			  ("resources,l", value<vector <string> >(&resources)->composing(), "Resource list")
		  			  ("name,N", value<string>(&name)->default_value("DefaultName"), "Job name")
		  			  ("variables,v", value<vector <string> >(&variables)->composing(), "Variables list to export")
		  			  ("exportAll,V", value<bool>(&exportAll)->implicit_value(true)->default_value(false), "Export all variables from context")
		  			  ("dependencies,W", value<vector <string> >(&dependencies)->composing(), "Add dependencies between jobs")
		  			  ;

		// Hidden options, will be allowed both on command line and in config file, but will not be shown to the user.
		options_description hidden("Hidden options");
		hidden.add_options()
					  ("input-file", value(&input_file),"input file")
					  ;

		options_description cmdline_options;
		cmdline_options.add(generic).add(config).add(hidden);

		options_description config_file_options;
		config_file_options.add(config).add(hidden);

		options_description visible("Allowed options");
		visible.add(generic).add(config);

		positional_options_description p;
		p.add("input-file", 1);

		variables_map vm;
		store(command_line_parser(ac, av).
				options(cmdline_options).positional(p).run(), vm);

		ifstream ifs("config.cfg");
		store(parse_config_file(ifs, config_file_options), vm);
		notify(vm);

		/*
		 * Checking options
		 */

		if (vm.count("help")) {
			cerr << visible << endl << "Mandatory option: " << endl << "<script>" << endl << endl;
			return EXIT_SUCCESS;
		}

		if (vm.count("version")) {
			cerr << "Queue Interpreter, version 1.0\n";
			return EXIT_SUCCESS;
		}

		if (vm.count("resources") && verbose)
		{
			cerr << "resources are: " << resources << "\n";
		}

		if (vm.count("input-file") && verbose)
		{
			cerr << "Input file is: " << input_file << endl;
		}

		map <string, string> maptron;
		if (vm.count("variables"))
		{
			if ( verbose ) cerr << "variables are: " << flush;
			for(vector<string>::iterator it = variables.begin(); it != variables.end(); ++it) {
				string text = *it;
				boost::char_separator<char> sep(",");
				boost::tokenizer< boost::char_separator<char> > tokens(text, sep);
				BOOST_FOREACH ( const string& temp_str, tokens) {
					int equal_pos=temp_str.find('=');
					string id = temp_str.substr(0,equal_pos);
					string val= temp_str.substr(equal_pos+1);
					if ( maptron[id] != "" ) { cerr << id << " defined twice" << endl; return EXIT_FAILURE; }
					maptron[id] = val;
					if ( verbose ) cerr << id << "=" << val << " " << flush;
				}
			}
			if ( verbose ) cerr << endl;
		}
		map<string,string>::const_iterator mit (maptron.begin()), mend(maptron.end());
		for(;mit!=mend;++mit) {
			setenv(mit->first.c_str(), mit->second.c_str(), true);
		}

		/*
		 * Running the script
		 */
		if ( !input_file.empty() ) {
			// Executing the input file
			string cmd = string("bash ").append(input_file).append(" 2>&1");
			if ( verbose ) cerr << "cmd to execute: " << endl << cmd << endl;

			FILE *pipe;
			pipe = popen ( cmd.c_str(), "r" );
			string data;
		    const int max_buffer = 256;
		    char buffer[max_buffer];
		    if (pipe) {
		        while (!feof(pipe))
		            if (fgets(buffer, max_buffer, pipe) != NULL)
		            	data.append(buffer);
		        returnVal = pclose(pipe);
		    }

		    if ( verbose ) cerr << data << endl;
		    ofstream out(name.append(".o").append(time_milli).c_str());
		    out << data;
		    out.close();
		}
		else {
			cout << "Undefined input file" << endl << endl;
			cout << visible << endl << "Mandatory option: " << endl << "<script>" << endl << endl;
			return EXIT_FAILURE;
		}
	}
	catch(exception& e)
	{
		cerr << e.what() << endl;
		return EXIT_FAILURE;
	}
	return returnVal;
}

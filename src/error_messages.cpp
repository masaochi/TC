// [namespace error_messages]
// print error messages and terminate the program

#include "include/header.hpp"

// Note! error messages are printed in std::cout

void error_messages::compare_list(const std::string &keyword, const std::string &input, const std::vector<std::string> &list)
{
    for (auto s : list)  // check whether "input" is included in "list"
    {
        if (input==s) return;  // found! do nothing
    }
    inappropriate_argument(keyword, input, "");  // no match then stop the prgoram with the error message
}    

// "condition": condition that "keyword" should satisfy    
void error_messages::inappropriate_argument(const std::string &keyword, const std::string &input, const std::string &condition)
{
    std::string s = "Error: inappropriate argument for " + keyword + "\nYour input: " + input;
    if (condition!="") { s += "\nCondition: " + condition; }
    stop(s);
}
void error_messages::inappropriate_argument(const std::string &keyword, const bool &input, const std::string &condition)
{
    std::string input_string = input ? "true" : "false";
    std::string s = "Error: inappropriate argument for " + keyword + "\nYour input: " + input_string;
    if (condition!="") { s += "\nCondition: " + condition; }
    stop(s);
}
void error_messages::inappropriate_argument(const std::string &keyword, const int &input, const std::string &condition)
{
    //    std::string s = "Error: inappropriate name for " + keyword + "\nYour input: " + std::to_string(input);
    std::ostringstream oss;
    oss << input;
    std::string s = "Error: inappropriate argument for " + keyword + "\nYour input: " + oss.str();
    if (condition!="") { s += "\nCondition: " + condition; }
    stop(s);
}
void error_messages::inappropriate_argument(const std::string &keyword, const double &input, const std::string &condition)
{
    //    std::string s = "Error: inappropriate name for " + keyword + "\nYour input: " + std::to_string(input);
    std::ostringstream oss;
    oss << input;
    std::string s = "Error: inappropriate argument for " + keyword + "\nYour input: " + oss.str();
    if (condition!="") { s += "\nCondition: " + condition; }
    stop(s);
}

void error_messages::not_found(const std::string &what, const std::string &where)
{
    stop("Error: could not find " + what + " in " + where);
}

void error_messages::cannot_open(const std::string &file_name)
{
    stop("Error: cannot open " + file_name);
}

void error_messages::stop(const std::string &message)
{
    std::cout << message << std::endl;
    MPI_Abort(MPI_COMM_WORLD,1);
}


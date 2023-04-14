// [namespace error_messages]
// print error messages and terminate the program

#ifndef TC_ERROR_MESSAGES_HPP
#define TC_ERROR_MESSAGES_HPP

namespace error_messages
{
    
// Call "inappropriate_argument()" when "input" is not included in "list"    
void compare_list(const std::string &keyword, const std::string &input, const std::vector<std::string> &list);

// Print an error message: an inappropriate argument is specified for "keyword", and stop the program
//   optional: "condition" is also displayed if condition!=""
void inappropriate_argument(const std::string &keyword, const std::string &input, const std::string &condition);
void inappropriate_argument(const std::string &keyword, const int &input, const std::string &condition);
void inappropriate_argument(const std::string &keyword, const double &input, const std::string &condition);

// Print an error message: "what" (keyword) is not found in "where" (file name), and stop the program
void not_found(const std::string &what, const std::string &where);

// Print an error message: cannot open "file_name", and stop the program
void cannot_open(const std::string &file_name);

// Stop the program with printing an error message "message"
void stop(const std::string &message);

}

#endif // TC_ERROR_MESSAGES_HPP


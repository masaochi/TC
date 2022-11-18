// [class MyClock]
// time (clock) information

#include "include/header.hpp"

MyClock::MyClock()
{
    start_ = std::chrono::system_clock::now();
    save_ = start_;
}

void MyClock::print_total_time(std::ostream *ost) const
{
    auto now = std::chrono::system_clock::now();
    std::string str = "Total time";
    print_time(ost, start_, now, str);
}

void MyClock::print_time_from_save(std::ostream *ost, const std::string name)
{
    auto now = std::chrono::system_clock::now();
    const std::string str = "Time for " + name;
    print_time(ost, save_, now, str);
    save_ = now;
}

void MyClock::print_time(std::ostream *ost, 
                         const std::chrono::system_clock::time_point &from,
                         const std::chrono::system_clock::time_point &to,
                         const std::string &str) const
{
    auto sec = std::chrono::duration_cast<std::chrono::seconds>(to - from).count();
    *ost << " " << str << " = " << (sec/86400) << " days " << ((sec%86400)/3600) << " h " << ((sec%3600)/60) << " m " <<  (sec%60) << " s" << std::endl;
    //*ost << std::endl;
}

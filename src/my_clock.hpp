// [class MyClock]
// time (clock) information

#ifndef TC_MY_CLOCK_HPP
#define TC_MY_CLOCK_HPP

class MyClock
{
private:
    std::chrono::system_clock::time_point start_;
    std::chrono::system_clock::time_point save_; // previous time point
    void print_time(std::ostream *ost,
                    const std::chrono::system_clock::time_point &from, 
                    const std::chrono::system_clock::time_point &to,
                    const std::string &str) const;

public:
    MyClock();
    void print_total_time(std::ostream *ost) const;
    void print_time_from_save(std::ostream *ost, const std::string name);
};

#endif // TC_MY_CLOCK_HPP

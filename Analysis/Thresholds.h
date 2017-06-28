#ifndef THRESHOLDS_H
#define THRESHOLDS_H

struct Thresholds_t
{
    // In  nanoseconds
    const double reject_time;
    const double delxtalk_reject_time;
    
    // In  % of pe
    const double after_pulse; 
    const double direct_xtalk;
    const double xtalk;
    const double time_dist;

} Thresholds {4, 2, 0.5, 1.17, 0.85, 0.4};

#endif
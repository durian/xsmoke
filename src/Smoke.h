#ifndef _SMOKE_H
#define _SMOKE_H

#include <string>

// ----------------------------------------------------------------------------
// Class
// ----------------------------------------------------------------------------

namespace SMOKE {

  class Smoke {

  public:
    std::string          prefsfilename;

    int   max_bols      = 50000; // if we have BOL
    float refresh = -1.0; // if we have BOL

    int   emit_rate = 120;
    float emit_velocity = 2.0;
    float flow = 0.90;
    int   throttle = 0;
    float movement = 0.40;
    float rgb_left[3];
    float rgb_right[3];
    float draw_delay = 0.0;
    float transparency_from = 0.4;
    float transparency_to = 0.01;
    float size_from = 0.0;
    float size_to = 4.0;
    float size_grow = 1.0;
    float emitter_from = 4.0;
    float emitter_to   = 5.0;
    float emitter_step = 2.0;
    // 
    float linger_min  = 30.0;
    float linger_time = 60.0; // random factor
    float linger_grow =  0.2;
    float colour_var  =  0.1;
    //
    float x_offset = 0.0; // left/right
    float y_offset = 0.0; // under or above
    float z_offset = 0.0; // behind
    //
    float start_vy  = 0.1; // start rise/fall
    float target_vy = 0.1; // target rise/fall
    //
    float wind_factor = 1.0;

    Smoke() {};

    ~Smoke() {};
  };
  
  //extern Global G;
}
#endif

#ifndef _GLOBAL_H
#define _GLOBAL_H

#include <string>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <fstream>
#include <vector>
#include <map>

#include "XPLMDataAccess.h"

#include "Smoke.h"

// ----------------------------------------------------------------------------
// Class
// ----------------------------------------------------------------------------

namespace SMOKE {

  class FloatWindow;

  std::string &ltrim(std::string &);
  std::string &rtrim(std::string &);
  std::string &trim(std::string &);
  void listify(const std::string&, std::vector<std::string>&);
    
  class Global {

  public:
    std::string          prefsfilename;
    std::string          keyfilename;
    bool                 loaded;
    bool                 valid_config;
        
    double elapsed;
    float dialog_interval_time;

    XPLMKeyFlags  gFlags = 0;
    XPLMKeyFlags  gPrevFlags = 0;
    unsigned char gVirtualKey = 0;
    unsigned char gChar = 0;

    int   max_bols      = 50000; // if we have BOL
    float smoke_refresh = -1.0; // if we have BOL

    int   plane_smoke_emit_rate = 120;
    float plane_smoke_emit_velocity = 2.0;
    float plane_smoke_flow = 0.90;
    int   plane_smoke_throttle = 0;
    float plane_smoke_movement = 0.40;
    float plane_smoke_rgb_left[3];
    float plane_smoke_rgb_right[3];
    float plane_smoke_draw_delay = 0.0;
    float plane_smoke_transparency_from = 0.4;
    float plane_smoke_transparency_to = 0.01;
    float plane_smoke_size_from = 0.0;
    float plane_smoke_size_to = 4.0;
    float plane_smoke_size_grow = 1.0;
    float plane_smoke_emitter_from = 4.0;
    float plane_smoke_emitter_to   = 5.0;
    float plane_smoke_emitter_step = 2.0;
    // 
    float plane_smoke_linger_min  = 30.0;
    float plane_smoke_linger_time = 60.0; // random factor
    float plane_smoke_linger_grow =  0.2;
    float plane_smoke_colour_var  =  0.1;
    //
    float plane_smoke_x_offset = 0.0; // left/right
    float plane_smoke_y_offset = 0.0; // under or above
    float plane_smoke_z_offset = 0.0; // behind
    //
    float plane_smoke_start_vy  = 0.1; // start rise/fall
    float plane_smoke_target_vy = 0.1; // target rise/fall
    //
    float plane_smoke_wind_factor = 1.0;
    //
    int scan_mode = 0; // 0 = loop, 1 = scan engine
    int ai_scan_mode = 0; // 0 = loop, 1 = scan engine
    int emit_num = 0;
    std::vector<float> emit_x;
    std::vector<float> emit_y;
    std::vector<float> emit_z;

    std::vector<double> mp_xs; // hold positions of multiplayer planes
    std::vector<double> mp_ys;
    std::vector<double> mp_zs;
    std::vector<float> mp_vxs; // hold positions of multiplayer planes
    std::vector<float> mp_vys;
    std::vector<float> mp_vzs;
    std::vector<float> mp_lat; // lat, lon
    std::vector<float> mp_lon;

    std::vector<float> cjs_ac_asl;
    std::vector<float> cjs_speed_kias;
    std::vector<float> cjs_heading_degT;
    std::vector<float> cjs_ac_lat;
    std::vector<float> cjs_ac_lon;

    float  nearest_ap_lat; // nearest airport lat, lon
    float  nearest_ap_lon;
    double nearest_ap_x; // nearest airport opengl
    double nearest_ap_y;
    double nearest_ap_z;

    std::string vkey;
    
    std::map<std::string, FloatWindow*> fwindows;

    // Guess what for
    Smoke ai_smoke;
    Smoke wt_smoke;
    
    // Constructor.
    Global();
    
    // Destructor.
    ~Global();

    bool read_prefs();
    bool read_prefs(const std::string&);
    bool read_smoke_prefs( const std::string& , Smoke& );
      
    float random_variation(int);
    float random_range(float,float);
    std::string random_string(size_t);

    float p_rand();
    float pn_rand();
    void interpolate_wind(float, float&, float&);
    void init_mp_datarefs();
    int get_mp_positions();
    void get_nearest_ap(double, double, float&, float&);
    int get_cjs();
  };
  
  extern Global G;
}
#endif

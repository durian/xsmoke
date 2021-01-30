#include <string>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <fstream>
#include <string>
#include <ctime>
#include <algorithm> 
#include <functional> 
#include <cctype>
#include <locale>
#include <random>

#include <stdlib.h>

#include "XPLMUtilities.h"
#include "XPLMDataAccess.h"
#include "XPLMPlanes.h"
#include "XPLMNavigation.h"
#include "XPLMScenery.h"
#include "XPLMGraphics.h"

#include "Global.h"
#include "Smoke.h"
#include "Log.h"
#include "dataref.h"

namespace SMOKE {

  DataRef<float> dr_wind_altitude_msl_m0("sim/weather/wind_altitude_msl_m[0]");
  DataRef<float> dr_wind_direction_degt0("sim/weather/wind_direction_degt[0]");
  DataRef<float> dr_wind_speed_kt0("sim/weather/wind_speed_kt[0]");
  DataRef<float> dr_wind_altitude_msl_m1("sim/weather/wind_altitude_msl_m[1]");
  DataRef<float> dr_wind_direction_degt1("sim/weather/wind_direction_degt[1]");
  DataRef<float> dr_wind_speed_kt1("sim/weather/wind_speed_kt[1]");
  DataRef<float> dr_wind_altitude_msl_m2("sim/weather/wind_altitude_msl_m[2]");
  DataRef<float> dr_wind_direction_degt2("sim/weather/wind_direction_degt[2]");
  DataRef<float> dr_wind_speed_kt2("sim/weather/wind_speed_kt[2]");

  XPLMDataRef cjs_num_ac = XPLMFindDataRef("cjs/world_traffic/num_aircraft");
  //DataRef<std::vector<float>> cjs_ac_asl("cjs/world_traffic/alt_asl"); // Altitudes of all aircraft (float array)
  //DataRef<std::vector<float>> cjs_speed_kias("cjs/world_traffic/speed_kias"); // Speeds of all aircraft (float array)
  //DataRef<std::vector<float>> cjs_heading_degT("cjs/world_traffic/heading_degT"); // Heading of all aircraft 
  //DataRef<std::vector<float>> cjs_ac_lat("cjs/world_traffic/aircraft_lat"); // Position of all aircraft (float array)
  //DataRef<std::vector<float>> cjs_ac_lon("cjs/world_traffic/aircraft_lon");
  XPLMDataRef cjs_ac_asl_dr = XPLMFindDataRef("cjs/world_traffic/alt_asl");
  XPLMDataRef cjs_speed_kias_dr = XPLMFindDataRef("cjs/world_traffic/speed_kias"); // Speeds of all aircraft
  XPLMDataRef cjs_heading_degT_dr = XPLMFindDataRef("cjs/world_traffic/heading_degT"); // Heading of all aircraft 
  XPLMDataRef cjs_ac_lat_dr = XPLMFindDataRef("cjs/world_traffic/aircraft_lat"); // Position of all aircraft
  XPLMDataRef cjs_ac_lon_dr = XPLMFindDataRef("cjs/world_traffic/aircraft_lon");
  
  //gPlaneX = XPLMGetDatad(gPlaneXDataRef);
  
  std::vector<XPLMDataRef> dr_mp_xs; // double
  std::vector<XPLMDataRef> dr_mp_ys;
  std::vector<XPLMDataRef> dr_mp_zs;
  std::vector<XPLMDataRef> dr_mp_vxs; // float
  std::vector<XPLMDataRef> dr_mp_vys;
  std::vector<XPLMDataRef> dr_mp_vzs;
  std::vector<XPLMDataRef> dr_mp_lat;
  std::vector<XPLMDataRef> dr_mp_lon;

  // ----------------------------------------------------------------------------
  // Code
  // ----------------------------------------------------------------------------

  // trim from start
  std::string &ltrim(std::string &s) {
    s.erase(s.begin(), std::find_if(s.begin(), s.end(), std::not1(std::ptr_fun<int, int>(std::isspace))));
    return s;
  }
  
  // trim from end
  std::string &rtrim(std::string &s) {
    s.erase(std::find_if(s.rbegin(), s.rend(), std::not1(std::ptr_fun<int, int>(std::isspace))).base(), s.end());
    return s;
  }
  
  // trim from both ends
  std::string &trim(std::string &s) {
    return ltrim(rtrim(s));
  }
  
  int parse_int(const std::string& s) { 
    int n;
    std::istringstream(s) >> n;
    return n;
  }

  void listify(const std::string& s, std::vector<std::string>& v) {
    char delim = ',';
    std::stringstream ss(s);
    std::string item;
    while (std::getline(ss, item, delim)) {
      trim(item);
      v.push_back(item);
    }
  }

  // Default values if no config file found.
  Global::Global() {
    plane_smoke_rgb_left[0] = 88;
    plane_smoke_rgb_left[1] = 0;
    plane_smoke_rgb_left[2] = 0;
    
    plane_smoke_rgb_right[0] = 0;
    plane_smoke_rgb_right[1] = 0;
    plane_smoke_rgb_right[2] = 88;
  }

  Global::~Global() {
    lg.xplm( "Deleted Global.\n" );
  }

  // Rereads the prefs if prefsfilename is defined
  bool Global::read_prefs() {
    lg.xplm( "read_prefs("+prefsfilename+")\n") ;
    if ( prefsfilename != "" ) {
      return read_prefs( prefsfilename );
    }
    return false;
  }
  
  bool Global::read_prefs( const std::string& filename ) {

    std::ifstream file( filename.c_str() );
    if ( ! file ) {
      lg.xplm( "WARNING: can't read prefs file. ("+filename+")\n" );
      return false;
    }

    bool keep_reading = true;

    scan_mode = 1; // scan, unless we find an emitter statement here
    // Some defaults
    G.plane_smoke_throttle = 0;
    
    std::string a_line;
    while( keep_reading ) {
      if ( ! std::getline( file, a_line ) ) {
	keep_reading = false;
	break;
      }
      if ( a_line.length() == 0 ) {
	continue;
      }
      if ( a_line.at(0) == '#' ) {
	continue;
      }
	 
      size_t pos = a_line.find( '=', 0 );
      if ( pos != std::string::npos ) {
	std::string lhs = a_line.substr( 0, pos );
	trim(lhs);
	std::string rhs = a_line.substr( pos+1 );
	trim(rhs);
	if ( (lhs != "") && (rhs != "") ) {
	  std::string tmp = lhs +"="+rhs+"\n";
	  lg.xplm( tmp );
	  // Some globals
	  if ( lhs == "smoke_max" ) {
	    int sm = int(std::stoi(rhs));
	    if ( (sm <= 0) || (sm > 500000) ) {
	      sm = 10000;
	    }
	    G.max_bols = sm;
	    lg.xplm("Smoke bols = "+std::to_string(G.max_bols)+"\n");
	  }
	  if ( lhs == "smoke_refresh" ) {
	    float sr = std::stof(rhs);
	    G.smoke_refresh = sr;
	    lg.xplm("Smoke refresh = "+std::to_string(G.smoke_refresh)+"\n");
	  }
	  if ( lhs == "plane_smoke_emit_rate" ) {
	    int ser = std::stoi(rhs);
	    if ( (ser <= 0) || (ser > 5000) ) {
	      ser = 120; 
	    }
	    G.plane_smoke_emit_rate = ser;
	    lg.xplm("Plane smoke emit rate: "+std::to_string(G.plane_smoke_emit_rate)+"\n" );
	  } else if ( lhs == "plane_smoke_emit_velocity" ) {
	    int sev = std::stoi(rhs);
	    if ( (sev < 0) || (sev > 100) ) {
	      sev = 2.0; 
	    }
	    G.plane_smoke_emit_velocity = sev;
	    lg.xplm("Plane smoke emit velocity: "+std::to_string(G.plane_smoke_emit_velocity)+"\n" );
	  } else if ( lhs == "plane_smoke_throttle" ) { // needs to be reset if changing to one w/o this param
	    int pst = std::stoi(rhs);
	    if ( pst != 1 ) {
	      pst = 0; 
	    }
	    G.plane_smoke_throttle = pst;
	    lg.xplm("Plane smoke throttle: "+std::to_string(G.plane_smoke_throttle)+"\n" );
	  } else if ( lhs == "plane_smoke_flow" ) {
	    float sf = std::stof(rhs);
	    if ( (sf <= 0) || (sf > 1.00) ) {
	      sf = 0.90; 
	    }
	    G.plane_smoke_flow = sf; // smoke flow of 0.99 is fastest; rand() > 0.01
	    lg.xplm("Plane smoke flow: "+std::to_string(G.plane_smoke_flow)+"\n" );
	  } else if ( lhs == "plane_smoke_movement" ) {
	    float sm = fabs(std::stof(rhs));
	    if ( sm > 100.0 ) {
	      sm = 100.00; 
	    }
	    G.plane_smoke_movement = sm; 
	    lg.xplm("Plane smoke movement: "+std::to_string(G.plane_smoke_movement)+"\n" );
	  } else if ( lhs == "plane_smoke_wind_factor" ) {
	    float wf = std::stof(rhs);
	    G.plane_smoke_wind_factor = wf; 
	    lg.xplm("Plane smoke wind factor: "+std::to_string(G.plane_smoke_wind_factor)+"\n" );
	  } else if ( lhs == "plane_smoke_colour_variation" ) {
	    float cv = std::stof(rhs);
	    if ( (cv < 0) || (cv > 1.0) ) {
	      cv = 0.1; 
	    }
	    G.plane_smoke_colour_var = cv; 
	    lg.xplm("Plane smoke colour variation: "+std::to_string(G.plane_smoke_colour_var)+"\n" );
	  } else if ( lhs == "plane_smoke_emit_offset" ) {
	    std::vector<std::string> bits;
	    listify( rhs, bits );
	    if ( bits.size() >= 3 ) {
	      float ox = 0.0; // x is left/right
	      float oy = 0.0; // y is under
	      float oz = 0.0; // z is behind
	      try {
		ox = std::stof(bits[0]);
		oy = std::stof(bits[1]);
		oz = std::stof(bits[2]);
	      } catch (...) {
		lg.xplm( "Can't read offset values, using 0, 0, 0 defaults.\n" );
	      }
	      G.plane_smoke_x_offset = ox;
	      G.plane_smoke_y_offset = oy;
	      G.plane_smoke_z_offset = oz;
	      lg.xplm("Plane smoke emit offset: "+std::to_string(G.plane_smoke_x_offset)+","+std::to_string(G.plane_smoke_y_offset)+","+std::to_string(G.plane_smoke_z_offset)+"\n" );
	    }
	    // if plane_smoke_emitter is not present, we also need scan_mode to be 1!
	  } else if ( lhs == "plane_smoke_emitter" ) { 
	    std::vector<std::string> bits;
	    listify( rhs, bits );
	    if ( bits.size() >= 3 ) {
	      float ef = 4.0;
	      float et = 5.0;
	      float es = 2.0;
	      try {
		ef = fabs(std::stof(bits[0]));
		et = fabs(std::stof(bits[1]));
		es = fabs(std::stof(bits[2]));
	      } catch (...) {
		lg.xplm( "Can't read emitter values, using 4, 5, 2 defaults.\n" );
	      }
	      if ( es <= 0.01 ) { // safety
		es = 1.0;
	      }
	      if ( et <= ef ) { // safety
		et = ef + 0.1;
	      }
	      G.plane_smoke_emitter_from = ef;
	      G.plane_smoke_emitter_to   = et;
	      G.plane_smoke_emitter_step = es;
	      scan_mode = 0;
	      lg.xplm("Plane smoke emitter: "+std::to_string(G.plane_smoke_emitter_from)+","+std::to_string(G.plane_smoke_emitter_to)+","+std::to_string(G.plane_smoke_emitter_step)+"\n" );
	    } else { // not 3 arguments
	      scan_mode = 1;  // scan datarefs instead
	      lg.xplm( "Will scan engine location datarefs.\n" );
	    }
	  } else if ( lhs == "plane_smoke_linger" ) { // time and growth
	    std::vector<std::string> bits;
	    listify( rhs, bits );
	    if ( bits.size() >= 3 ) {
	      float lm = 30.0;
	      float lt = 60.0;
	      float lr =  1.0;
	      try {
		lm = fabs(std::stof(bits[0]));
		lt = fabs(std::stof(bits[1]));
		lr = std::stof(bits[2]); // can be negative
	      } catch (...) {
		lg.xplm( "Can't read linger values, using 30, 60, 1 defaults.\n" );
	      }
	      G.plane_smoke_linger_min  = lm;
	      G.plane_smoke_linger_time = lt;
	      G.plane_smoke_linger_grow = lr;
	      lg.xplm("Plane smoke linger: "+std::to_string(G.plane_smoke_linger_min)+","+std::to_string(G.plane_smoke_linger_time)+","+std::to_string(G.plane_smoke_linger_grow)+"\n" );
	    }
	  } else if ( lhs == "plane_smoke_rgb_left" ) { // left and right!
	    std::vector<std::string> bits;
	    listify( rhs, bits );
	    if ( bits.size() >= 3 ) {
	      float r = 0.0;
	      float g = 0.0;
	      float b = 0.0;
	      try {
		r = fabs(std::stof(bits[0]));
		g = fabs(std::stof(bits[1]));
		b = fabs(std::stof(bits[2]));
	      } catch (...) {
		lg.xplm( "Can't read rgb values, using 0,0,0 defaults.\n" );
	      }
	      G.plane_smoke_rgb_left[0] = r/255.0;
	      G.plane_smoke_rgb_left[1] = g/255.0;
	      G.plane_smoke_rgb_left[2] = b/255.0;
	      lg.xplm("Plane smoke RGB left: "+std::to_string(G.plane_smoke_rgb_left[0])+","+std::to_string(G.plane_smoke_rgb_left[1])+","+std::to_string(G.plane_smoke_rgb_left[2])+"\n" );
	    }
	  } else if ( lhs == "plane_smoke_rgb_right" ) { // left and right!
	    std::vector<std::string> bits;
	    listify( rhs, bits );
	    if ( bits.size() >= 3 ) {
	      float r = 0.0;
	      float g = 0.0;
	      float b = 0.0;
	      try {
		r = fabs(std::stof(bits[0]));
		g = fabs(std::stof(bits[1]));
		b = fabs(std::stof(bits[2]));
	      } catch (...) {
		lg.xplm( "Can't read rgb values, using 0,0,0 defaults.\n" );
	      }
	      G.plane_smoke_rgb_right[0] = r/255.0;
	      G.plane_smoke_rgb_right[1] = g/255.0;
	      G.plane_smoke_rgb_right[2] = b/255.0;
	      lg.xplm("Plane smoke RGB right: "+std::to_string(G.plane_smoke_rgb_right[0])+","+std::to_string(G.plane_smoke_rgb_right[1])+","+std::to_string(G.plane_smoke_rgb_right[2])+"\n" );
	    }
	  } else if ( lhs == "plane_smoke_draw_delay" ) {
	    float dd = fabs(std::stof(rhs));
	    if ( dd > 600.0 ) {
	      dd = 600.0;  //do we need a limit?
	    }
	    G.plane_smoke_draw_delay = dd;
	    lg.xplm("Plane smoke draw delay: "+std::to_string(G.plane_smoke_draw_delay)+"\n" );
	  } else if ( lhs == "plane_smoke_target_vy" ) {
	    float ty = std::stof(rhs);
	    G.plane_smoke_target_vy = ty;
	    lg.xplm("Plane smoke target vy: "+std::to_string(G.plane_smoke_target_vy)+"\n" );
	  } else if ( lhs == "plane_smoke_vy" ) {
	    std::vector<std::string> bits;
	    listify( rhs, bits );
	    if ( bits.size() >= 2 ) {
	      float sy = 1.0;
	      float ty = 1.0;
	      try {
		sy = std::stof(bits[0]);
		ty = std::stof(bits[1]);
	      } catch (...) {
		lg.xplm( "Can't read vy values, using 1, 1 defaults.\n" );
	      }
	      G.plane_smoke_start_vy  = sy;
	      G.plane_smoke_target_vy = ty;
	      lg.xplm("Plane smoke vy: "+std::to_string(G.plane_smoke_start_vy)+","+std::to_string(G.plane_smoke_target_vy)+"\n" );
	    }
	  }  else if ( lhs == "plane_smoke_transparency" ) {
	    std::vector<std::string> bits;
	    listify( rhs, bits );
	    if ( bits.size() >= 2 ) {
	      float tr_from = 0.2;
	      float tr_to = 0.0;
	      try {
		tr_from = fabs(std::stof(bits[0]));
		tr_to   = fabs(std::stof(bits[1]));
	      } catch (...) {
		lg.xplm( "Can't read transparency values, using 0.2,0.0 defaults.\n" );
	      }
	      G.plane_smoke_transparency_from = tr_from;
	      G.plane_smoke_transparency_to   = tr_to;
	      lg.xplm("Plane smoke transparency: "+std::to_string(G.plane_smoke_transparency_from)+","+std::to_string(G.plane_smoke_transparency_to)+"\n" );
	    }
	  } else if ( lhs == "plane_smoke_size" ) {
	    std::vector<std::string> bits;
	    listify( rhs, bits );
	    if ( bits.size() >= 3 ) {
	      float si_from = 0.02;
	      float si_to   = 0.20;
	      float si_grow = 0.02; //m/s, diameter
	      try {
		si_from = fabs(std::stof(bits[0]));
		si_to   = fabs(std::stof(bits[1]));
		si_grow = fabs(std::stof(bits[2]));
		if ( si_grow <= 0.001 ) {
		  lg.xplm("Positive growth needed.\n");
		  si_grow = 0.001;
		}
		if ( si_to < si_from ) {
		  lg.xplm("Second value must be larger than first.\n");
		  si_to = si_from + 1.0;
		}
	      } catch (...) {
		lg.xplm( "Can't read size values, using 0.02,0.2,0.02 defaults.\n" );
	      }
	      G.plane_smoke_size_from = si_from;
	      G.plane_smoke_size_to   = si_to;
	      G.plane_smoke_size_grow = si_grow;
	      lg.xplm("Plane smoke size: "+std::to_string(G.plane_smoke_size_from)+","+std::to_string(G.plane_smoke_size_to)+","+std::to_string(G.plane_smoke_size_grow)+"\n" );
	    }
	  }
	}
      }
    }

    
    file.close();
    prefsfilename = filename;

    // if all invalid, fill by hand with flare/return false?
    return true;
  }

  // This should be the only one FIXME
  bool Global::read_smoke_prefs( const std::string& filename, Smoke& S ) {

    std::ifstream file( filename.c_str() );
    if ( ! file ) {
      lg.xplm( "WARNING: can't read prefs file. ("+filename+")\n" );
      return false;
    }

    bool keep_reading = true;

    ai_scan_mode = 1; // scan, unless we find an emitter statement here PJB FIXME THIS CHANGES USERPLANE SMOKE!
    // Some defaults
    S.throttle = 0;
    
    std::string a_line;
    while( keep_reading ) {
      if ( ! std::getline( file, a_line ) ) {
	keep_reading = false;
	break;
      }
      if ( a_line.length() == 0 ) {
	continue;
      }
      if ( a_line.at(0) == '#' ) {
	continue;
      }
	 
      size_t pos = a_line.find( '=', 0 );
      if ( pos != std::string::npos ) {
	std::string lhs = a_line.substr( 0, pos );
	trim(lhs);
	std::string rhs = a_line.substr( pos+1 );
	trim(rhs);
	if ( (lhs != "") && (rhs != "") ) {
	  std::string tmp = lhs +"="+rhs+"\n";
	  lg.xplm( tmp );
	  // Some globals
	  if ( lhs == "smoke_max" ) {
	    int sm = int(std::stoi(rhs));
	    if ( (sm <= 0) || (sm > 500000) ) {
	      sm = 10000;
	    }
	    S.max_bols = sm;
	    lg.xplm("Smoke bols = "+std::to_string(S.max_bols)+"\n");
	  }
	  if ( lhs == "smoke_refresh" ) {
	    float sr = std::stof(rhs);
	    S.refresh = sr;
	    lg.xplm("Smoke refresh = "+std::to_string(S.refresh)+"\n");
	  }
	  if ( lhs == "plane_smoke_emit_rate" ) {
	    int ser = std::stoi(rhs);
	    if ( (ser <= 0) || (ser > 5000) ) {
	      ser = 120; 
	    }
	    S.emit_rate = ser;
	    lg.xplm("Plane smoke emit rate: "+std::to_string(S.emit_rate)+"\n" );
	  } else if ( lhs == "plane_smoke_emit_velocity" ) {
	    int sev = std::stoi(rhs);
	    if ( (sev < 0) || (sev > 100) ) {
	      sev = 2.0; 
	    }
	    S.emit_velocity = sev;
	    lg.xplm("Plane smoke emit velocity: "+std::to_string(S.emit_velocity)+"\n" );
	  } else if ( lhs == "plane_smoke_throttle" ) { // needs to be reset if changing to one w/o this param
	    int pst = std::stoi(rhs);
	    if ( pst != 1 ) {
	      pst = 0; 
	    }
	    S.throttle = pst;
	    lg.xplm("Plane smoke throttle: "+std::to_string(S.throttle)+"\n" );
	  } else if ( lhs == "plane_smoke_flow" ) {
	    float sf = std::stof(rhs);
	    if ( (sf <= 0) || (sf > 1.00) ) {
	      sf = 0.90; 
	    }
	    S.flow = sf; // smoke flow of 0.99 is fastest; rand() > 0.01
	    lg.xplm("Plane smoke flow: "+std::to_string(S.flow)+"\n" );
	  } else if ( lhs == "plane_smoke_movement" ) {
	    float sm = fabs(std::stof(rhs));
	    if ( sm > 100.0 ) {
	      sm = 100.00; 
	    }
	    S.movement = sm; 
	    lg.xplm("Plane smoke movement: "+std::to_string(S.movement)+"\n" );
	  } else if ( lhs == "plane_smoke_wind_factor" ) {
	    float wf = std::stof(rhs);
	    S.wind_factor = wf; 
	    lg.xplm("Plane smoke wind factor: "+std::to_string(S.wind_factor)+"\n" );
	  } else if ( lhs == "plane_smoke_colour_variation" ) {
	    float cv = std::stof(rhs);
	    if ( (cv < 0) || (cv > 1.0) ) {
	      cv = 0.1; 
	    }
	    S.colour_var = cv; 
	    lg.xplm("Plane smoke colour variation: "+std::to_string(S.colour_var)+"\n" );
	  } else if ( lhs == "plane_smoke_emit_offset" ) {
	    std::vector<std::string> bits;
	    listify( rhs, bits );
	    if ( bits.size() >= 3 ) {
	      float ox = 0.0; // x is left/right
	      float oy = 0.0; // y is under
	      float oz = 0.0; // z is behind
	      try {
		ox = std::stof(bits[0]);
		oy = std::stof(bits[1]);
		oz = std::stof(bits[2]);
	      } catch (...) {
		lg.xplm( "Can't read offset values, using 0, 0, 0 defaults.\n" );
	      }
	      S.x_offset = ox;
	      S.y_offset = oy;
	      S.z_offset = oz;
	      lg.xplm("Plane smoke emit offset: "+std::to_string(S.x_offset)+","+std::to_string(S.y_offset)+","+std::to_string(S.z_offset)+"\n" );
	    }
	    // if plane_smoke_emitter is not present, we also need scan_mode to be 1!
	  } else if ( lhs == "plane_smoke_emitter" ) { 
	    std::vector<std::string> bits;
	    listify( rhs, bits );
	    if ( bits.size() >= 3 ) {
	      float ef = 4.0;
	      float et = 5.0;
	      float es = 2.0;
	      try {
		ef = fabs(std::stof(bits[0]));
		et = fabs(std::stof(bits[1]));
		es = fabs(std::stof(bits[2]));
	      } catch (...) {
		lg.xplm( "Can't read emitter values, using 4, 5, 2 defaults.\n" );
	      }
	      if ( es <= 0.01 ) { // safety
		es = 1.0;
	      }
	      if ( et <= ef ) { // safety
		et = ef + 0.1;
	      }
	      S.emitter_from = ef;
	      S.emitter_to   = et;
	      S.emitter_step = es;
	      ai_scan_mode = 0;
	      lg.xplm("Plane smoke emitter: "+std::to_string(S.emitter_from)+","+std::to_string(S.emitter_to)+","+std::to_string(S.emitter_step)+"\n" );
	    } else { // not 3 arguments
	      ai_scan_mode = 1;  // scan datarefs instead
	      lg.xplm( "Will scan engine location datarefs.\n" );
	    }
	  } else if ( lhs == "plane_smoke_linger" ) { // time and growth
	    std::vector<std::string> bits;
	    listify( rhs, bits );
	    if ( bits.size() >= 3 ) {
	      float lm = 30.0;
	      float lt = 60.0;
	      float lr =  1.0;
	      try {
		lm = fabs(std::stof(bits[0]));
		lt = fabs(std::stof(bits[1]));
		lr = std::stof(bits[2]); // can be negative
	      } catch (...) {
		lg.xplm( "Can't read linger values, using 30, 60, 1 defaults.\n" );
	      }
	      S.linger_min  = lm;
	      S.linger_time = lt;
	      S.linger_grow = lr;
	      lg.xplm("Plane smoke linger: "+std::to_string(S.linger_min)+","+std::to_string(S.linger_time)+","+std::to_string(S.linger_grow)+"\n" );
	    }
	  } else if ( lhs == "plane_smoke_rgb_left" ) { // left and right!
	    std::vector<std::string> bits;
	    listify( rhs, bits );
	    if ( bits.size() >= 3 ) {
	      float r = 0.0;
	      float g = 0.0;
	      float b = 0.0;
	      try {
		r = fabs(std::stof(bits[0]));
		g = fabs(std::stof(bits[1]));
		b = fabs(std::stof(bits[2]));
	      } catch (...) {
		lg.xplm( "Can't read rgb values, using 0,0,0 defaults.\n" );
	      }
	      S.rgb_left[0] = r/255.0;
	      S.rgb_left[1] = g/255.0;
	      S.rgb_left[2] = b/255.0;
	      lg.xplm("Plane smoke RGB left: "+std::to_string(S.rgb_left[0])+","+std::to_string(S.rgb_left[1])+","+std::to_string(S.rgb_left[2])+"\n" );
	    }
	  } else if ( lhs == "plane_smoke_rgb_right" ) { // left and right!
	    std::vector<std::string> bits;
	    listify( rhs, bits );
	    if ( bits.size() >= 3 ) {
	      float r = 0.0;
	      float g = 0.0;
	      float b = 0.0;
	      try {
		r = fabs(std::stof(bits[0]));
		g = fabs(std::stof(bits[1]));
		b = fabs(std::stof(bits[2]));
	      } catch (...) {
		lg.xplm( "Can't read rgb values, using 0,0,0 defaults.\n" );
	      }
	      S.rgb_right[0] = r/255.0;
	      S.rgb_right[1] = g/255.0;
	      S.rgb_right[2] = b/255.0;
	      lg.xplm("Plane smoke RGB right: "+std::to_string(S.rgb_right[0])+","+std::to_string(S.rgb_right[1])+","+std::to_string(S.rgb_right[2])+"\n" );
	    }
	  } else if ( lhs == "plane_smoke_draw_delay" ) {
	    float dd = fabs(std::stof(rhs));
	    if ( dd > 600.0 ) {
	      dd = 600.0;  //do we need a limit?
	    }
	    S.draw_delay = dd;
	    lg.xplm("Plane smoke draw delay: "+std::to_string(S.draw_delay)+"\n" );
	  } else if ( lhs == "plane_smoke_target_vy" ) {
	    float ty = std::stof(rhs);
	    S.target_vy = ty;
	    lg.xplm("Plane smoke target vy: "+std::to_string(S.target_vy)+"\n" );
	  } else if ( lhs == "plane_smoke_vy" ) {
	    std::vector<std::string> bits;
	    listify( rhs, bits );
	    if ( bits.size() >= 2 ) {
	      float sy = 1.0;
	      float ty = 1.0;
	      try {
		sy = std::stof(bits[0]);
		ty = std::stof(bits[1]);
	      } catch (...) {
		lg.xplm( "Can't read vy values, using 1, 1 defaults.\n" );
	      }
	      S.start_vy  = sy;
	      S.target_vy = ty;
	      lg.xplm("Plane smoke vy: "+std::to_string(S.start_vy)+","+std::to_string(S.target_vy)+"\n" );
	    }
	  }  else if ( lhs == "plane_smoke_transparency" ) {
	    std::vector<std::string> bits;
	    listify( rhs, bits );
	    if ( bits.size() >= 2 ) {
	      float tr_from = 0.2;
	      float tr_to = 0.0;
	      try {
		tr_from = fabs(std::stof(bits[0]));
		tr_to   = fabs(std::stof(bits[1]));
	      } catch (...) {
		lg.xplm( "Can't read transparency values, using 0.2,0.0 defaults.\n" );
	      }
	      S.transparency_from = tr_from;
	      S.transparency_to   = tr_to;
	      lg.xplm("Plane smoke transparency: "+std::to_string(S.transparency_from)+","+std::to_string(S.transparency_to)+"\n" );
	    }
	  } else if ( lhs == "plane_smoke_size" ) {
	    std::vector<std::string> bits;
	    listify( rhs, bits );
	    if ( bits.size() >= 3 ) {
	      float si_from = 0.02;
	      float si_to   = 0.20;
	      float si_grow = 0.02; //m/s, diameter
	      try {
		si_from = fabs(std::stof(bits[0]));
		si_to   = fabs(std::stof(bits[1]));
		si_grow = fabs(std::stof(bits[2]));
		if ( si_grow <= 0.001 ) {
		  lg.xplm("Positive growth needed.\n");
		  si_grow = 0.001;
		}
		if ( si_to < si_from ) {
		  lg.xplm("Second value must be larger than first.\n");
		  si_to = si_from + 1.0;
		}
	      } catch (...) {
		lg.xplm( "Can't read size values, using 0.02,0.2,0.02 defaults.\n" );
	      }
	      S.size_from = si_from;
	      S.size_to   = si_to;
	      S.size_grow = si_grow;
	      lg.xplm("Plane smoke size: "+std::to_string(S.size_from)+","+std::to_string(S.size_to)+","+std::to_string(S.size_grow)+"\n" );
	    }
	  }
	}
      }
    }

    
    file.close();
    prefsfilename = filename;

    // if all invalid, fill by hand with flare/return false?
    return true;
  }

  float Global::random_variation(int variation) {
    // return 0.9-1.1 for random_variation(10) -> a multiplication factor
    float x = float((rand() % (2*variation) + (100-variation)) / (float)100.0);
    return x;
  }
  
  float Global::random_range(float l, float h) {
    static std::random_device rd;     // only used once to initialise (seed) engine
    std::mt19937 rng(rd());    // random-number engine used (Mersenne-Twister in this case)
    std::uniform_real_distribution<float> dis(l,h); // guaranteed unbiased
    auto r = dis(rng);
    return r;
  }

  std::string Global::random_string( size_t length ) {
    auto randchar = []() -> char {
      const char charset[] =
      "0123456789"
      "ABCDEFGHIJKLMNOPQRSTUVWXYZ"
      "abcdefghijklmnopqrstuvwxyz";
      const size_t max_index = (sizeof(charset) - 1);
      return charset[ rand() % max_index ];
    };
    std::string str(length,0);
    std::generate_n( str.begin(), length, randchar );
    return str;
  }
  
  float Global::p_rand() { // positive random 0..1
    return float(std::rand()) / float(RAND_MAX);
  }

  float Global::pn_rand() { // positive-negative random -1..1
    return (float(std::rand())-float(std::rand())) / float(RAND_MAX);
  }

  void Global::interpolate_wind(float msl, float& vel, float& dir) {
    // why this happens?
    if ( ( (dr_wind_altitude_msl_m1 - dr_wind_altitude_msl_m0) <= 0 )
	 ||
	 ( (dr_wind_altitude_msl_m2 - dr_wind_altitude_msl_m1) <= 0 )
	 ) {
      //lg.xplm( "Wind levels at strange altitudes\n" ); // this spams too much
      vel = 0;
      dir = 0;
      return;
    }	     
    if ( msl <= dr_wind_altitude_msl_m0 ) {
      // return layer 0
      vel = dr_wind_speed_kt0;
      dir = dr_wind_direction_degt0;
    } else if ( msl <= dr_wind_altitude_msl_m1 ) {
      // interpolate between 0-1
      float factor =
	(msl - dr_wind_altitude_msl_m0)
	/
	(dr_wind_altitude_msl_m1 - dr_wind_altitude_msl_m0);
      //lg.xplm( "FACTOR "+std::to_string(factor)+"\n" );
      vel = dr_wind_speed_kt0 + ( factor * ( dr_wind_speed_kt1 - dr_wind_speed_kt0 ) );
      dir = dr_wind_direction_degt0 + ( factor * ( dr_wind_direction_degt1 - dr_wind_direction_degt0 ) );
    } else {
      float factor =
	(msl - dr_wind_altitude_msl_m1)
	/
	(dr_wind_altitude_msl_m2 - dr_wind_altitude_msl_m1);
      vel = dr_wind_speed_kt1 + ( factor * ( dr_wind_speed_kt2 - dr_wind_speed_kt1 ) );
      dir = dr_wind_direction_degt1 + ( factor * ( dr_wind_direction_degt2 - dr_wind_direction_degt1 ) );
    }

    vel *= 0.51444; // convert knots to m/s

    // take layer 0
    vel = dr_wind_speed_kt0 * 0.51444; // m/s
    dir = dr_wind_direction_degt0;
  }

  void Global::init_mp_datarefs() {
    char buffer [128];
    int n;
    dr_mp_xs.clear();
    dr_mp_ys.clear();
    dr_mp_zs.clear();
    dr_mp_vxs.clear();
    dr_mp_vys.clear();
    dr_mp_vzs.clear();
    dr_mp_lat.clear();
    dr_mp_lon.clear();

    for( int i = 1; i < 20; i++ ) { // datarefs for all
      n = sprintf(buffer, "sim/multiplayer/position/plane%d_x", i);
      lg.xplm( std::string(buffer)+"\n" );
      //dr_mp_xs.push_back( XPLMFindDataRef("sim/multiplayer/position/plane1_x") );
      XPLMDataRef x = XPLMFindDataRef(buffer);
      if ( x ) {
	dr_mp_xs.push_back( x );
      } else {
	lg.xplm( "ERROR reading dr "+std::to_string(i)+"\n" );
      }
      n = sprintf(buffer, "sim/multiplayer/position/plane%d_v_x", i);
      lg.xplm( std::string(buffer)+"\n" );
      x = XPLMFindDataRef(buffer);
      if ( x ) {
	dr_mp_vxs.push_back( x );
      } else {
	lg.xplm( "ERROR reading dr "+std::to_string(i)+"\n" );
      }
      //
      n = sprintf(buffer, "sim/multiplayer/position/plane%d_y", i);
      x = XPLMFindDataRef(buffer);
      if ( x ) {
	dr_mp_ys.push_back( x );
      } else {
	lg.xplm( "ERROR reading dr "+std::to_string(i)+"\n" );
      }
      n = sprintf(buffer, "sim/multiplayer/position/plane%d_v_y", i);
      x = XPLMFindDataRef(buffer);
      if ( x ) {
	dr_mp_vys.push_back( x );
      } else {
	lg.xplm( "ERROR reading dr "+std::to_string(i)+"\n" );
      }
      //
      n = sprintf(buffer, "sim/multiplayer/position/plane%d_z", i);
      x = XPLMFindDataRef(buffer);
      if ( x ) {
	dr_mp_zs.push_back( x );
      } else {
	lg.xplm( "ERROR reading dr "+std::to_string(i)+"\n" );
      }
      n = sprintf(buffer, "sim/multiplayer/position/plane%d_v_z", i);
      x = XPLMFindDataRef(buffer);
      if ( x ) {
	dr_mp_vzs.push_back( x );
      } else {
	lg.xplm( "ERROR reading dr "+std::to_string(i)+"\n" );
      }
      //
      //
      n = sprintf(buffer, "sim/multiplayer/position/plane%d_lat", i);
      x = XPLMFindDataRef(buffer);
      if ( x ) {
	dr_mp_lat.push_back( x );
      } else {
	lg.xplm( "ERROR reading dr "+std::to_string(i)+"\n" );
      }
      n = sprintf(buffer, "sim/multiplayer/position/plane%d_lon", i);
      x = XPLMFindDataRef(buffer);
      if ( x ) {
	dr_mp_lon.push_back( x );
      } else {
	lg.xplm( "ERROR reading dr "+std::to_string(i)+"\n" );
      }
    }
    // NOTE THAT PLANE1 IS AT POSITION 0 !
    lg.xplm( "POSITIONS: "+std::to_string(dr_mp_xs.size())+"\n" );
  }

  //void Global::get_mp_positions( std::vector<double>& xs, std::vector<double>& ys, std::vector<double>& zs ) {
  int Global::get_mp_positions() {
    int	mp_planes;
    XPLMCountAircraft( &mp_planes, 0, 0 );
    --mp_planes;
    //lg.xplm( "mp_planes: "+std::to_string(mp_planes)+"\n" ); // AI planes 1, mp_planes = 2
    
    G.mp_xs.clear();
    G.mp_ys.clear();
    G.mp_zs.clear();
    G.mp_vxs.clear();
    G.mp_vys.clear();
    G.mp_vzs.clear();
    G.mp_lat.clear();
    G.mp_lon.clear();
    
    for( int i = 0; i < mp_planes; i++ ) { // trying just the first 5
      //lg.xplm( "Getting "+std::to_string(i)+"\n" );
      double mp_x = XPLMGetDatad( dr_mp_xs[i] );
      G.mp_xs.push_back( mp_x );
      float mp_vx = XPLMGetDataf( dr_mp_vxs[i] );
      G.mp_vxs.push_back( mp_vx );
      
      double mp_y = XPLMGetDatad( dr_mp_ys[i] );
      G.mp_ys.push_back( mp_y );
      float mp_vy = XPLMGetDataf( dr_mp_vys[i] );
      G.mp_vys.push_back( mp_vy );
      
      double mp_z = XPLMGetDatad( dr_mp_zs[i] );
      G.mp_zs.push_back( mp_z );
      float mp_vz = XPLMGetDataf( dr_mp_vzs[i] );
      G.mp_vzs.push_back( mp_vz );
      
      double mp_lat = XPLMGetDatad( dr_mp_lat[i] );
      G.mp_lat.push_back( mp_lat );
      double mp_lon = XPLMGetDatad( dr_mp_lon[i] );
      G.mp_lon.push_back( mp_lon );

      //std::string m =" POS AI1 "+std::to_string(ai_x)+", "+std::to_string(ai_y)+", "+std::to_string(ai_z)+"\n";
      //lg.xplm(m);
    }
    /*
    int	planeCount;
    XPLMCountAircraft( &planeCount, 0, 0 );
    char FileName[256], AircraftPath[256];
    char* pAircraft[256];
    for (long index = 0; index < planeCount; ++index) {
      XPLMGetNthAircraftModel(index, FileName, AircraftPath);
      pAircraft[index] = (char *)AircraftPath;
      lg.xplm( std::string((char*)AircraftPath)+"\n");
      if (XPLMAcquirePlanes((char **)&pAircraft, NULL, NULL)) {
	//XPLMSetAircraftModel(index, AircraftPath);
	//m =" POS AI1 "+std::to_string(dr_mp_lx)+", "+std::to_string(dr_mp_ly)+", "+std::to_string(dr_mp_lz)+"\n";
	//XPLMDisableAIForPlane(index);
	//std::string m = " POS AI "+std::to_string(dr_mp_lx)+", "+std::to_string(dr_mp_ly)+", "+std::to_string(dr_mp_lz)+"\n";
	//lg.xplm(m);
      }
    }
    */
    return mp_planes;
  }

  void Global::get_nearest_ap(double plane_lat, double plane_lon, float& latitude, float& longitude) {
    float lat = static_cast<float>(plane_lat);
    float lon = static_cast<float>(plane_lon);
    
    XPLMNavRef closest_ap = XPLMFindNavAid(
					   NULL, //const char *         inNameFragment,    /* Can be NULL */
					   NULL, //const char *         inIDFragment,    /* Can be NULL */
					   (float*)&lat, 
					   (float*)&lon, //float *              inLon,    /* Can be NULL */
					   NULL, //int *                inFrequency,    /* Can be NULL */
					   xplm_Nav_Airport);
    
    char id[32];
    char name[256];
    XPLMGetNavAidInfo(closest_ap, NULL, &latitude, &longitude, NULL, NULL, NULL, id, name, NULL);
    lg.xplm( "NEAREST Plane lat, lon: "+std::to_string(lat)+", "+std::to_string(lon)+"\n");
    lg.xplm( "NEAREST AP: "+std::string(id)+", "+std::string(name)+", "+
	     std::to_string(latitude)+", "+std::to_string(longitude)+"\n");

    G.nearest_ap_lat = latitude;
    G.nearest_ap_lon = longitude;

    XPLMWorldToLocal( G.nearest_ap_lat, G.nearest_ap_lon, 0,
		      &G.nearest_ap_x, &G.nearest_ap_y, &G.nearest_ap_z );

    XPLMProbeRef hProbe;
    hProbe = XPLMCreateProbe(xplm_ProbeY);
    XPLMProbeInfo_t info = { 0 };
    info.structSize = sizeof(info);
    // If we have a hit then return Y coordinate
    if ( XPLMProbeTerrainXYZ( hProbe, G.nearest_ap_x, G.nearest_ap_y, G.nearest_ap_z, &info) == xplm_ProbeHitTerrain ) {
      lg.xplm( "G.nearest_ap_y="+std::to_string(G.nearest_ap_y)+", info.locationY="+std::to_string(info.locationY)+"\n" );
      G.nearest_ap_y = info.locationY;
    }
  }

  // Fills a global vector...
  int Global::get_cjs() {
    if ( ! cjs_num_ac ) { // No WT3 plugin
      return 0;
    }
    
    int num_ac = XPLMGetDatai(cjs_num_ac);
    
    //lg.xplm( "Number of cjs ac: "+std::to_string(num_ac)+"\n" );

    cjs_ac_asl.reserve(num_ac);
    cjs_ac_asl = {};
    cjs_speed_kias.reserve(num_ac);
    cjs_speed_kias = {};
    cjs_heading_degT.reserve(num_ac);
    cjs_heading_degT = {};
    cjs_ac_lat.reserve(num_ac);
    cjs_ac_lat = {};
    cjs_ac_lon.reserve(num_ac);
    cjs_ac_lon = {};

    for( int i = 0; i < num_ac; i++ ) {
      cjs_ac_lat[i] = 0;
      cjs_ac_lon[i] = 0;
      cjs_ac_asl[i] = 0;
      cjs_speed_kias[i] = 0;
    }
    
    XPLMGetDatavf(cjs_ac_asl_dr, &cjs_ac_asl[0], 0, num_ac);
    XPLMGetDatavf(cjs_speed_kias_dr, &cjs_speed_kias[0], 0, num_ac);
    XPLMGetDatavf(cjs_heading_degT_dr, &cjs_heading_degT[0], 0, num_ac);
    XPLMGetDatavf(cjs_ac_lat_dr, &cjs_ac_lat[0], 0, num_ac);
    XPLMGetDatavf(cjs_ac_lon_dr, &cjs_ac_lon[0], 0, num_ac);

    /*
    for( int i = 0; i < num_ac; i++ ) {
      float lat = cjs_ac_lat[i];
      float lon = cjs_ac_lon[i];
      float asl = cjs_ac_asl[i];
      float spd = cjs_speed_kias[i];
      lg.xplm("CJS: "+std::to_string(i)+", "+std::to_string(lat)+", "+std::to_string(lon)+", "+std::to_string(asl)+", "+std::to_string(spd)+"\n");
    }
    */
    /*
    for( int i = 0; i < num_ac; i++ ) {
      lg.xplm("1\n");
      float lat = cjs_ac_lat[i];
      lg.xplm("2\n");
      float lon = cjs_ac_lon[i];
      lg.xplm("3\n");
      float asl = cjs_ac_asl[i];
      lg.xplm("4\n");
      lg.xplm("Engine: "+std::to_string(lat)+", "+std::to_string(lon)+", "+std::to_string(asl)+"\n");
    }
    */

    return num_ac;
  }

  Global G;
}

// The End --------

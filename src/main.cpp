/*
  (c) pjb 2015, 2016, 2017, 2018, 2019, 2020, 2021
  -----BEGIN GEEK CODE BLOCK-----
  Version: 3.1
  GCS/H/IT/P d--?pu s+:+ a++ C++$ UL++$ P++ L+++$ E++$ W++$ N+++ !o-- K++ !w-- O M+ V PS++ PE- Y+ PGP t++ 5++ X++ R tv+ b++ DI+ D++ G e++++ h r++ y+
  ------END GEEK CODE BLOCK------ 
  https://www.joereiss.net/geek/geek2.cgi

21:24/21:25 wrong screen

make lin
cp ../src/build-lin/liblin.xpl ~/xplane11/Resources/plugins/xsmoke/64/lin.xpl

*/
#if IBM
#include <windows.h>
BOOL APIENTRY DllMain( HANDLE hModule,
                       DWORD  ul_reason_for_call,
                       LPVOID lpReserved
                     )
{
    switch (ul_reason_for_call)
    {
        case DLL_PROCESS_ATTACH:
        case DLL_THREAD_ATTACH:
        case DLL_THREAD_DETACH:
        case DLL_PROCESS_DETACH:
            break;
    }
    return TRUE;
}
#endif

#if APL
#include <OpenGL/OpenGL.h>
#include <OpenGL/glu.h>
#elif IBM
#define WIN32_LEAN_AND_MEAN
#include <windows.h>
//#include "gl/glew.h"
#include <GL/gl.h>
#include <GL/glu.h>
#elif LIN
//#include <GL/glew.h>
#include <GL/gl.h>
#include <GL/glu.h>
#endif

#define TRUE 1
#define FALSE 0

#define M_PI           3.14159265358979323846

#include <string>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <fstream>
#include <algorithm>
#include <chrono>
#include <vector>

#include <stdio.h>
#include <string.h>
#include <time.h>
#include <math.h>

#include "XPLMProcessing.h"
#include "XPLMDataAccess.h"
#include "XPLMUtilities.h"
#include "XPLMPlugin.h"
#include "XPLMMenus.h"
#include "XPLMGraphics.h"
#include "XPLMScenery.h"
#include "XPLMCamera.h"
#include "XPLMPlanes.h"
#include "XPLMNavigation.h"
#include "XPLMDisplay.h"

#include "dataref.h"
#include "main.h"
#include "Global.h"
#include "Log.h"
#include "Bol.h"

/*
#define PERLIN 1
#include "PerlinNoise.h"
#include "ppm.h"
*/

using namespace SMOKE;

// 0.97 adds ground test on downwards only
//           smoke 180 degrees from plane psi
// 0.98 smoke throttle, plane smoke.ini
std::string VERSION = "1.1.2";

DataRef<int>    dr_sim_paused("sim/time/paused");
DataRef<int>    dr_sim_speed("sim/time/sim_speed");
DataRef<float>  dr_sim_trts("sim/time/total_running_time_sec");
DataRef<float>  dr_sim_tfts("sim/time/total_flight_time_sec");
DataRef<int>    dr_is_reverse_float_z("sim/graphics/view/is_reverse_float_z");
  
static bool smokemode   = true;
static bool aismokemode = true;
static bool wtsmokemode = true;
static int  headmode    = -1; // 0 or larger is index to closest mp
static int  findapmode  = -1; // 0 or larger is index to closest mp

DataRef<double> dr_plane_lx("sim/flightmodel/position/local_x");
DataRef<double> dr_plane_ly("sim/flightmodel/position/local_y");
DataRef<double> dr_plane_lz("sim/flightmodel/position/local_z");

DataRef<float> dr_plane_lvx("sim/flightmodel/position/local_vx");
DataRef<float> dr_plane_lvy("sim/flightmodel/position/local_vy");
DataRef<float> dr_plane_lvz("sim/flightmodel/position/local_vz");

DataRef<double> dr_plane_lon("sim/flightmodel/position/longitude");
DataRef<double> dr_plane_lat("sim/flightmodel/position/latitude");

/*
sim/weather/wind_now_x_msc	float	n	meters/sec	Wind direction vector in OpenGL coordinates, X component.
sim/weather/wind_now_y_msc	float	n	meters/sec	Wind direction vector in OpenGL coordinates, Y component.
sim/weather/wind_now_z_msc	float	n	meters/sec	Wind direction vector in OpenGL coordinates, Z component.
*/
DataRef<float> dr_wind_now_x_msc("sim/weather/wind_now_x_msc");
DataRef<float> dr_wind_now_y_msc("sim/weather/wind_now_y_msc");
DataRef<float> dr_wind_now_z_msc("sim/weather/wind_now_z_msc");
/*
wind_direction_degt	float	660+	no	[0-359)	The effective direction of the wind at the plane's location.
wind_speed_kt	float	660+	no	kts >= 0	The effective speed of the wind at the plane's location.
*/
DataRef<float> dr_wind_direction_degt("sim/weather/wind_direction_degt");
DataRef<float> dr_wind_speed_kt("sim/weather/wind_speed_kt");

DataRef<float>  dr_true_airspeed("sim/flightmodel/position/true_airspeed"); // m/s
DataRef<double> dr_plane_elevation("sim/flightmodel/position/elevation");
DataRef<float>  dr_plane_y_agl("sim/flightmodel/position/y_agl");
float reference_h = 0.0;

DataRef<float> dr_plane_ogl_psi("sim/flightmodel/position/psi"); // psi? true_psi? mag_psi?
DataRef<float> dr_plane_psi("sim/flightmodel/position/true_psi"); // psi? true_psi? mag_psi?
DataRef<float> dr_plane_the("sim/flightmodel/position/true_theta");
DataRef<float> dr_plane_phi("sim/flightmodel/position/phi");

DataRef<std::vector<float>> dr_plane_q("sim/flightmodel/position/q");

DataRef<float> dr_ph_psi("sim/graphics/view/pilots_head_psi", ReadWrite); // Position of pilot's head heading
DataRef<float> dr_ph_the("sim/graphics/view/pilots_head_the", ReadWrite); // Position of pilot's head pitch
DataRef<float> dr_ph_phi("sim/graphics/view/pilots_head_phi", ReadWrite); // Position of the pilot's head roll
DataRef<float> dr_ph_x("sim/graphics/view/pilots_head_x", ReadWrite); // Position of pilot's head relative to CG, X
DataRef<float> dr_ph_y("sim/graphics/view/pilots_head_y", ReadWrite); // Position of pilot's head relative to CG, Y
DataRef<float> dr_ph_z("sim/graphics/view/pilots_head_z", ReadWrite); // Position of pilot's head relative to CG, Z

// sim/flightmodel2/engines/
//Datarefs for the physical engines, props and rotors.
//location_x_mtr	float[engines]	900+	no	meters	Engine location, meters x, y, z, with respect to the default center of gravity.
//location_y_mtr	float[engines]	900+	no	meters 	Engine location, meters x, y, z, with respect to the default center of gravity.
//location_z_mtr	float[engines]	900+	no	meters 	Engine location, meters x, y, z, with respect to the default center of gravity.
bool rescan = true;
DataRef<int> dr_acf_num_engines("sim/aircraft/engine/acf_num_engines");
DataRef<std::vector<float>> dr_location_x_mtr("sim/flightmodel2/engines/location_x_mtr");
DataRef<std::vector<float>> dr_location_y_mtr("sim/flightmodel2/engines/location_y_mtr");
DataRef<std::vector<float>> dr_location_z_mtr("sim/flightmodel2/engines/location_z_mtr");

float last_tfts = 0.0; // to determine menus/maps

// Increase rate with throttle
DataRef<float> dr_throttle_ratio_all("sim/cockpit2/engine/actuators/throttle_ratio_all");

double plane_x_prev;
double plane_y_prev;
double plane_z_prev;

std::vector<CBol*>::iterator xcb_i; // for use in loop/erase

XPLMMenuID	myMenu;
int		mySubMenuItem;

enum gpxlog_status {MENU_TOGGLE, MENU_AI_TOGGLE, MENU_WT_TOGGLE, MENU_RELOAD};

XPLMCommandRef  SmokeCommand;
int SmokeCommandHandler(XPLMCommandRef, XPLMCommandPhase, void *);
XPLMCommandRef  AISmokeCommand;
int AISmokeCommandHandler(XPLMCommandRef, XPLMCommandPhase, void *);
XPLMCommandRef  WTSmokeCommand;
int WTSmokeCommandHandler(XPLMCommandRef, XPLMCommandPhase, void *);
XPLMCommandRef  EulerCommand;
int EulerCommandHandler(XPLMCommandRef, XPLMCommandPhase, void *);
XPLMCommandRef  FindAPCommand;
int FindAPCommandHandler(XPLMCommandRef, XPLMCommandPhase, void *);

float emit_subcounter    = 0.0;
int   max_particles_seen = 0;

// quaternion library
//http://forums.x-plane.org/index.php?/forums/topic/52920-plane-movement-and-physics-engine-toggling/

// too many globals
std::string the_dir;
std::vector<std::string> inifiles;
int current_inifile = 0; // actually, smoke.ini
int ai_inifile = 0; // pointer in vector
int wt_inifile = 0; // pointer in vector

#define DEG_TO_RAD_2 M_PI / 360.0
#ifndef M_PI
#define M_PI           3.14159265358979323846
#endif

std::string tail(std::string const& source, size_t const length) {
  if (length >= source.size()) { return source; }
  return source.substr(source.size() - length);
}

// distance between points, in meters.
static double distanceto(double lat0, double lon0, double lat1, double lon1) {
  double slat = sin((lat1-lat0) * (double)(M_PI/360));
  double slon = sin((lon1-lon0) * (double)(M_PI/360));
  double aa   = slat*slat + cos(lat0 * (double)(M_PI/180)) * cos(lat1 * (double)(M_PI/180)) * slon * slon;
  return 6378145.0 * 2 * atan2(sqrtf(aa), sqrt(1-aa));
}

int closest_mp() {
  double lat0 = dr_plane_lat;
  double lon0 = dr_plane_lon;
  int idx = -1;
  double min_d = 100000;
  int mp_planes = G.get_mp_positions();

  for( int i = 0; i < mp_planes; i++ ) {  
    double d = distanceto( lat0, lon0, G.mp_lat[i], G.mp_lon[i] );
    //lg.xplm( "DISTANCE: "+std::to_string(lat0)+", "+std::to_string(lon0)+", "+std::to_string(G.mp_lat[i])+", "+std::to_string(G.mp_lon[i])+"\n" );
    if ( d < min_d ) {
      idx = i;
      min_d = d;
    }
  }
  lg.xplm( "CLOSEST MP: "+std::to_string(idx)+", "+std::to_string(min_d)+"\n" );
  return idx;
}

// Convert plane coords / opengl coords, to get a point on e.g. engine in
// plane coords to world/opengl coords
void conversion(float x_plane, float y_plane, float z_plane, // = source location in airplane coordinates.  
		  float phi, float psi, float the, 
		  float local_x, float local_y, float local_z,  //plane's location in the world 
		  float& x_wrl, float& y_wrl, float& z_wrl) {
    // OUTPUTS:(x_wrl, y_wrl, z_wrl) = transformed location in world.
    float x_phi = x_plane*cos(phi) + y_plane*sin(phi);
    float y_phi = y_plane*cos(phi) - x_plane*sin(phi);
    float z_phi = z_plane;
    float x_the = x_phi;
    float y_the = y_phi*cos(the) - z_phi*sin(the);
    float z_the = z_phi*cos(the) + y_phi*sin(the);
    x_wrl = x_the*cos(psi) - z_the*sin(psi) + local_x;
    y_wrl = y_the                           + local_y;
    z_wrl = z_the*cos(psi) + x_the*sin(psi) + local_z;
  }

static XPLMProbeRef hProbe = XPLMCreateProbe(xplm_ProbeY);
float height(double x, double y, double z) {
  XPLMProbeInfo_t info = { 0 };
  info.structSize = sizeof(info);  
  // If we have a hit then return Y coordinate
  if (XPLMProbeTerrainXYZ( hProbe, x, y, z, &info) == xplm_ProbeHitTerrain) {
    return info.locationY;
  }
  return -1.0; // not smart, actual value
}

std::string rounded(float number) {
  std::stringstream ss;
  ss << std::fixed << std::setprecision(2) << number;
  return ss.str();
}

/*
sim/graphics/view/pilots_head_x float   y       meters  Position of pilot's head relative to CG, X
sim/graphics/view/pilots_head_y float   y       meters  Position of pilot's head relative to CG, Y
sim/graphics/view/pilots_head_z float   y       meters  Position of pilot's head relative to CG, Z

(2, -1, 3) and ( 1, 4, -3)

v = v0 - v1 = (1, -5, 6)
so: r = (2, -1, 3) + t(1, -5, 6 )
    x=2+t, y=-1+-5t, z=3+6t

*/
void euler(int idx) {
  int num_mp = G.get_mp_positions();
  if ( num_mp == 0 ) {
    return;
  }
  double ai_x = G.mp_xs[idx]; 
  double ai_y = G.mp_ys[idx];
  double ai_z = G.mp_zs[idx];
  //lg.xplm( "EULER AI: "+std::to_string(ai_x)+", "+std::to_string(ai_y)+", "+std::to_string(ai_z)+"\n" );

  double u_x = dr_plane_lx + dr_ph_x;
  double u_y = dr_plane_ly + dr_ph_y;
  double u_z = dr_plane_lz + dr_ph_z;
  //lg.xplm( "EULER US: "+std::to_string(u_x)+", "+std::to_string(u_y)+", "+std::to_string(u_z)+", "+std::to_string(dr_plane_psi)+"\n" );
  
  double delta_x = ai_x - u_x;
  double delta_y = ai_y - u_y;
  double delta_z = ai_z - u_z;

  // extrapolate, make new camera point (or head pos)
  // then aim in atan2 results, like head?
  /*
  double t     = 10.0; // meters in opengl NEEDS NORMALISATION!
  double l     = sqrt( pow(delta_x, 2) + pow(delta_y, 2) + pow(delta_z, 2) );
  double new_x = u_x + ( t * delta_x / l ); // positive is right/east
  double new_y = u_y + ( t * delta_y / l);
  double new_z = u_z + ( t * delta_z / l); // positive is back/south
  dr_ph_x = new_x; 
  dr_ph_y = new_y;
  dr_ph_z = new_z; // should only be done once...
  */
  
  // atan2, Returns the principal value of the arc tangent of y/x, expressed in radians.
  // atan2(double y, double x);
  // 0 degs is over the x-axis?
  
  float dir = atan2( delta_z, delta_x ); 
  dir = dir * 180.0 / 3.1415 + 90 - dr_plane_psi;
  if ( dir < 0 ) {
    dir = 360 + dir;
  }
  float pitch = atan2( delta_y, sqrt( pow(delta_x, 2) + pow(delta_z, 2) ) );
  pitch = pitch * 180.0 / 3.1415;
  //- dr_plane_the; // ADJUST dr_plane ?

  //lg.xplm( "EULER NU: "+std::to_string(dr_ph_psi)+", "+std::to_string(dr_ph_the)+", "+std::to_string(dr_ph_phi)+"\n" );
  //lg.xplm( "EULER NU: "+std::to_string(dr_ph_x)+", "+std::to_string(dr_ph_y)+", "+std::to_string(dr_ph_z)+"\n" );
  //lg.xplm( "EULER: "+std::to_string(dir)+", "+std::to_string(pitch)+"\n" );
  
  dr_ph_psi = dir;
  dr_ph_the = pitch;
  dr_ph_phi = -dr_plane_phi; //0; // or plane_phi adjustment ?
}

void move_head_towards(float lat, float lon, float ele) {
  double ai_x;
  double ai_y;
  double ai_z;
  XPLMWorldToLocal( lat, lon, ele, &ai_x, &ai_y, &ai_z ); // this can (is?) be cached
  
  double u_x = dr_plane_lx + dr_ph_x;
  double u_y = dr_plane_ly + dr_ph_y;
  double u_z = dr_plane_lz + dr_ph_z;
  
  double delta_x = ai_x - u_x;
  double delta_y = ai_y - u_y;
  double delta_z = ai_z - u_z;

  float dir = atan2( delta_z, delta_x ); 
  dir = dir * 180.0 / 3.1415 + 90 - dr_plane_psi;
  if ( dir < 0 ) {
    dir = 360 + dir;
  }
  float pitch = atan2( delta_y, sqrt( pow(delta_x, 2) + pow(delta_z, 2) ) );
  pitch = pitch * 180.0 / 3.1415;
  //pitch = pitch - dr_plane_the;  // PITCH NEEDS TO BE ANGLED, IF 90degs ANGLE, NONE

  dr_ph_psi = dir;
  dr_ph_the = pitch;
  dr_ph_phi = -dr_plane_phi; // or plane_phi adjustment ?
}

/*
  Init with engines, or loop settings
  Finds the location of the engines, where smoke will be emitted.
*/
void scan_engines() {

  lg.xplm("Scan mode: "+std::to_string(G.scan_mode)+"\n");
  if ( G.scan_mode == 1 ) {
    //float quat_0 = dr_plane_q[0];
    // dr_location_x_mtr
    int acf_num_engines = dr_acf_num_engines;
    
    G.emit_num = acf_num_engines;
    G.emit_x.clear();
    //G.emit_x.resize(acf_num_engines);
    G.emit_y.clear();
    //G.emit_y.resize(acf_num_engines);
    G.emit_z.clear();
    //G.emit_z.resize(acf_num_engines);
    
    lg.xplm("Num Engines (scan): "+rounded(G.emit_num)+"\n");
    for( int i = 0; i < acf_num_engines; i++ ) {
      float x = dr_location_x_mtr[i];
      float y = dr_location_y_mtr[i];
      float z = dr_location_z_mtr[i];
      lg.xplm("Engine: "+rounded(x)+", "+rounded(y)+", "+rounded(z)+"\n");
      /*
      G.emit_x[i] = x;
      G.emit_y[i] = y + G.plane_smoke_y_offset;
      G.emit_z[i] = z + G.plane_smoke_z_offset;
      */
      G.emit_x.push_back(x + G.plane_smoke_x_offset);
      G.emit_y.push_back(y + G.plane_smoke_y_offset);
      G.emit_z.push_back(z + G.plane_smoke_z_offset);
    }
  } else {
    G.emit_x.clear();
    //G.emit_x.resize(10);
    G.emit_y.clear();
    //G.emit_y.resize(10);
    G.emit_z.clear();
    //G.emit_z.resize(10);
    int c = 0;
    for ( float i = G.plane_smoke_emitter_from; i <= G.plane_smoke_emitter_to; i += G.plane_smoke_emitter_step ) {
      G.emit_x.push_back(i+G.plane_smoke_x_offset);
      G.emit_y.push_back(G.plane_smoke_y_offset);
      G.emit_z.push_back(G.plane_smoke_z_offset);
      lg.xplm("Engine: "+rounded(G.emit_x[c])+", "+rounded(G.emit_y[c])+", "+rounded(G.emit_z[c])+"\n");
      c += 1;
      if ( i > 0.0 ) {
	G.emit_x.push_back(-i-G.plane_smoke_x_offset); // note reversal
	G.emit_y.push_back(G.plane_smoke_y_offset);
	G.emit_z.push_back(G.plane_smoke_z_offset);
	lg.xplm("Engine: "+rounded(G.emit_x[c])+", "+rounded(G.emit_y[c])+", "+rounded(G.emit_z[c])+"\n");
	c += 1;
      }
    }
    G.emit_num = c;
    lg.xplm("Num Engines: "+rounded(G.emit_num)+"\n");
  }
}

// http://forums.x-plane.org/index.php?/forums/topic/49598-question-about-xplmgetdirectorycontents-function/
void scan_directory(std::string& the_dir, std::vector<std::string>& res) {
  char   FilesInFolder[15000];
  int    NumberOfFiles;
  int    TotalNumberOfFiles;
  char*  FileIndex[500];

  XPLMGetDirectoryContents(the_dir.c_str(), 0, FilesInFolder, sizeof(FilesInFolder), FileIndex, 500,
			   &TotalNumberOfFiles, &NumberOfFiles);
  lg.xplm( "Files: "+std::to_string(TotalNumberOfFiles)+", "+std::to_string(NumberOfFiles)+"\n" );
  std::string fn = "";
  int i = 0;
  char c;
  while ( NumberOfFiles > 0 ) {
    while ( (c = FilesInFolder[i]) != '\0' ) {
      fn += c;
      i += 1;
    }
    if ( (fn.size() > 4) && (fn.substr(fn.size()-4) == std::string(".ini")) ) { 
      lg.xplm( "File :"+fn+"\n" ); 
      res.push_back(fn);
    }
    fn = "";
    i += 1;
    NumberOfFiles--;
  }
}

// Find smoke.ini in aircraft directory
// ac/smoke.ini contains what... name to smoke? complete ini file?
int ac_directory(std::string& p) {
  int RetVal = 0;
  int Position;
  char FileName[512], AircraftPath[1024], SettingsFilePath[512];   
  // Get the current aircraft path.
  XPLMGetNthAircraftModel(0, FileName, AircraftPath);
  //lg.xplm( "AC0: "+std::string(AircraftPath)+"\n");
  //20150709 12:24:24.345: AC: /Volumes/Luna/X-Plane 10 clean/Aircraft/General Aviation/Columbia-400/c400.acf
  size_t found;
  p = std::string(AircraftPath);
  const char *sep = XPLMGetDirectorySeparator();
  found = p.find_last_of(sep, std::string::npos);
  if (found != std::string::npos) {
    p = p.substr(0,found)+sep;
    p += "smoke.ini";
    // if file exists, return 1
    std::ifstream file( p.c_str() );
    if ( ! file ) {
      return 0;
    }
    file.close();
    return 1;
  }
  return 0;
}
	
PLUGIN_API int XPluginStart(char *outName, char *outSig, char *outDesc) {
  strcpy(outName, "X-Smoke");
  strcpy(outSig,  "org.durian.xsmoke");
  strcpy(outDesc, "A plugin which generates smoke trails.");
  std::string compile_date("X-SMOKE plugin compiled on " __DATE__ " at " __TIME__ "\n");
  lg.xplm( compile_date.c_str() );
  lg.xplm( "Version "+VERSION+"\n" );

  int xp_ver;
  int xplm_ver;
  XPLMHostApplicationID xplm_host;
  XPLMGetVersions(&xp_ver, &xplm_ver, &xplm_host);
  lg.xplm( "xp_ver: "+std::to_string(xp_ver) + "\n" );
  
  time_t rawtime;
  struct tm * timeinfo;
  time ( &rawtime );
  // http://www.onlineconversion.com/unix_time.htm
  if ( rawtime >= (time_t)2546330088 ) { // Friday, September 9, 2050 9:54:48 AM
    lg.xplm("This demo version has expired\n");
    lg.xplm("Please update to a newer version.\n");
    return 0;
  }

  // ----------------------------------------------------------------------
#ifdef PERLIN
  /*
  perlin.precalc_perlins();
  for ( int i = 0; i < 100000; i++ ) {
    float n = perlin.next_perlin();
    if ( (n<0.2) || (n>0.8) ) {
      lg.xplm("PERLIN "+std::to_string(i)+"="+std::to_string(n)+"\n");
    }
  }
  */
  unsigned int width = 512, height = 512;
  
  // Create an empty PPM image
  ppm image(width, height);
  
  // Create a PerlinNoise object with the reference permutation vector
  PerlinNoise pn;
  
  unsigned int kk = 0;
  // Visit every pixel of the image and assign a color generated with Perlin noise
  for(unsigned int i = 0; i < height; ++i) {     // y
    for(unsigned int j = 0; j < width; ++j) {  // x
      double x = (double)j/((double)width);
      double y = (double)i/((double)height);
      
      // Typical Perlin noise
      //double n = pn.noise(8.0 * x, 8.0 * y, 0.0);
      double n = pn.octave_noise(x,y, 0.0, 4, 1); // octaves, persistence
      //lg.xplm("PERLIN "+std::to_string(x)+","+std::to_string(y)+"="+std::to_string(n)+"\n");
      
      // Wood like structure
      /*
	n = 20 * pn.noise(x, y, 0.8);
	n = n - floor(n);
      */
      // Map the values to the [0, 255] interval, for simplicity we use 
      // tones of grey
      image.r[kk] = floor(255 * n);
      image.g[kk] = floor(255 * n);
      image.b[kk] = floor(255 * n);
      kk++;
    }
  }
  // Save the image in a binary PPM file
  //image.write("figure_8_R.ppm");
#endif  
  // ----------------------------------------------------------------------


  // https://developer.x-plane.com/code-sample/instanced-drawing-sample/
  
  XPLMEnableFeature("XPLM_USE_NATIVE_PATHS", 1);  

  char filebase[255];
  XPLMGetSystemPath(filebase); // Locate the X-System directory
  const char *sep = XPLMGetDirectorySeparator();

  // The (default)initialisation file.
  /*
  std::string prefsfile = std::string(filebase) + "Resources" + sep + "plugins" + sep+ "smoke" + sep + "smoke.ini";  
  G.prefsfilename = prefsfile;
  G.read_prefs(G.prefsfilename);
  */

  // PJB FIXME crashes if not found:
  the_dir = std::string(filebase) + "Resources" + sep + "plugins" + sep+ "xsmoke";
  scan_directory(the_dir, inifiles);
    
  load_textures();
  
  // First we put a new menu item into the plugin menu
  mySubMenuItem = XPLMAppendMenuItem(XPLMFindPluginsMenu(), // Plugins menu 
				     "X-Smoke", // Menu title
				     0,  // Item ref
				     1); // English 

  // Now create a submenu attached to our menu items
  myMenu = XPLMCreateMenu("X-Smoke",
			  XPLMFindPluginsMenu(),
			  mySubMenuItem, /* Menu Item to attach to. */
			  MyMenuHandlerCallback,/* The handler */
			  0);			/* Handler Ref */


  XPLMAppendMenuItem(myMenu, "Toggle Smoke",  (void*)MENU_TOGGLE, 1);
  XPLMCheckMenuItem(myMenu, MENU_TOGGLE, xplm_Menu_Unchecked);

  XPLMAppendMenuItem(myMenu, "Toggle AI Smoke",  (void*)MENU_AI_TOGGLE, 1);
  XPLMCheckMenuItem(myMenu, MENU_AI_TOGGLE, xplm_Menu_Unchecked);

  XPLMAppendMenuItem(myMenu, "Toggle WT Smoke",  (void*)MENU_WT_TOGGLE, 1);
  XPLMCheckMenuItem(myMenu, MENU_WT_TOGGLE, xplm_Menu_Unchecked);

  XPLMAppendMenuItem(myMenu, "Rescan directory",  (void*)MENU_RELOAD, 1);

  XPLMAppendMenuSeparator(myMenu);
  
  // The above should be rescan-directory, and then a list with ini files
  //
  for( int i = 0; i < inifiles.size(); i++ ) {
    //for ( auto const &inifile : inifiles ) {
    //XPLMAppendMenuItem(myMenu, inifiles[i].c_str(),  (void*)((long)i+4), 1);
    XPLMAppendMenuItem(myMenu, inifiles[i].c_str(),  (void*)((uintptr_t)(i+5)), 1); // i offset
    if ( inifiles[i] == "smoke.ini" ) {
      current_inifile = i;
    }
    if ( inifiles[i] == "ai.ini" ) {
      ai_inifile = i;
    }
    if ( inifiles[i] == "wt.ini" ) {
      wt_inifile = i;
    }
    // PJB TODO read WT and AI.ini files here?
  }
  XPLMCheckMenuItem(myMenu, ((long)current_inifile+5), xplm_Menu_Checked);

  if ( inifiles[(long)current_inifile] != "smoke.ini" ) {
    G.read_prefs( the_dir + sep + inifiles[(long)current_inifile] );
  } else {
    std::string tmp;
    int r = ac_directory(tmp);
    if ( r == 1 ) {
      G.read_prefs( tmp );
    } else {
      G.read_prefs( the_dir + sep + inifiles[(long)current_inifile] );
    }
  }

  G.read_smoke_prefs( the_dir + sep + "ai.ini", G.ai_smoke );
  G.read_smoke_prefs( the_dir + sep + "wt.ini", G.wt_smoke );
  lg.xplm( "wt.ini rgb "+std::to_string(G.wt_smoke.rgb_right[0])+"\n" );
  
  XPLMRegisterFlightLoopCallback(DeferredInitNewAircraftFLCB, -1, NULL);

  SmokeCommand = XPLMCreateCommand("durian/xsmoke/toggle", "Toggle Smoke");
  XPLMRegisterCommandHandler(SmokeCommand,        // in Command name
			     SmokeCommandHandler, // in Handler
			     0,                  // Receive input before plugin windows. (or 1?)
			     (void *) 0);        // inRefcon.

  AISmokeCommand = XPLMCreateCommand("durian/xsmoke/toggleai", "Toggle AI Smoke");
  XPLMRegisterCommandHandler(AISmokeCommand,        // in Command name
			     AISmokeCommandHandler, // in Handler
			     0,                  // Receive input before plugin windows. (or 1?)
			     (void *) 0);        // inRefcon.

  WTSmokeCommand = XPLMCreateCommand("durian/xsmoke/togglewt", "Toggle WT Smoke");
  XPLMRegisterCommandHandler(WTSmokeCommand,        // in Command name
			     WTSmokeCommandHandler, // in Handler
			     0,                  // Receive input before plugin windows. (or 1?)
			     (void *) 0);        // inRefcon.

  EulerCommand = XPLMCreateCommand("durian/xsmoke/togglelook", "Look at closest AI plane");
  XPLMRegisterCommandHandler(EulerCommand,        // in Command name
			     EulerCommandHandler, // in Handler
			     0,                  // Receive input before plugin windows. (or 1?)
			     (void *) 0);        // inRefcon.

  FindAPCommand = XPLMCreateCommand("durian/xsmoke/findairport", "Look at closest airport");
  XPLMRegisterCommandHandler(FindAPCommand,        // in Command name
			     FindAPCommandHandler, // in Handler
			     0,                  // Receive input before plugin windows. (or 1?)
			     (void *) 0);        // inRefcon.

  ////XPLMRegisterDrawCallback(DrawBolBillboard, xplm_Phase_Modern3D, 0, NULL);
  XPLMRegisterDrawCallback(DrawBolBillboardVK, xplm_Phase_Modern3D, 0, NULL);

  
  //XPLMRegisterDrawCallback(DrawBolBillboard, xplm_Phase_Airplanes, 0, NULL);

  ////XPLMRegisterDrawCallback(DrawBolBillboard, xplm_Phase_Objects, 0, NULL);

  scan_engines();
  
  smokemode   = false;
  aismokemode = false;
  wtsmokemode = false;

  XPLMDataRef cjs_num_ac = XPLMFindDataRef("cjs/world_traffic/num_aircraft");
  if ( ! cjs_num_ac ) { // No WT3 plugin  
    XPLMEnableMenuItem(myMenu, MENU_WT_TOGGLE, 0);
    lg.xplm( "Disabled WT3 functionality.\n" );
  }

  if ( dr_is_reverse_float_z == 1 ) {
    lg.xplm( "Using Vulkan rendering.\n" );
  }
  
  lg.xplm( "Initialised.\n" );

  return 1;
}

PLUGIN_API void	XPluginStop(void) {
  smokemode = false;
  XPLMCheckMenuItem(myMenu, MENU_TOGGLE, xplm_Menu_Unchecked);
  aismokemode = false;
  XPLMCheckMenuItem(myMenu, MENU_AI_TOGGLE, xplm_Menu_Unchecked);
  wtsmokemode = false;
  XPLMCheckMenuItem(myMenu, MENU_WT_TOGGLE, xplm_Menu_Unchecked);

  ////XPLMUnregisterDrawCallback(DrawBolBillboard, xplm_Phase_Modern3D, 0, NULL);
  XPLMUnregisterDrawCallback(DrawBolBillboardVK, xplm_Phase_Modern3D, 0, NULL);
    
  //XPLMUnregisterDrawCallback(DrawBolBillboard, xplm_Phase_Airplanes, 0, NULL);

  XPLMUnregisterFlightLoopCallback(MyFlightLoopCallback, NULL);
  // remove menus?
}

PLUGIN_API int XPluginEnable(void) {
  lg.xplm("XPluginEnable\n");
  return 1;
}

PLUGIN_API void XPluginDisable(void) {
  smokemode = false;
  aismokemode = false;
  wtsmokemode = false;
  lg.xplm("XPluginDisable.\n");
}

PLUGIN_API void XPluginReceiveMessage(XPLMPluginID inFromWho, long inMessage, void *inParam ) {

  (void)inFromWho;
  (void)inParam;

  if ( inMessage == XPLM_MSG_PLANE_CRASHED ) { // 101
    smokemode = false;
    aismokemode = false;
    XPLMCheckMenuItem(myMenu, MENU_TOGGLE, xplm_Menu_Unchecked);
    XPLMCheckMenuItem(myMenu, MENU_AI_TOGGLE, xplm_Menu_Unchecked);
    rescan = true;
  }
  
  if ( inMessage == XPLM_MSG_PLANE_LOADED ) { // 102
    smokemode = false;
    aismokemode = false;
    XPLMCheckMenuItem(myMenu, MENU_TOGGLE, xplm_Menu_Unchecked);
    XPLMCheckMenuItem(myMenu, MENU_AI_TOGGLE, xplm_Menu_Unchecked);
    rescan = true;
    
    const char *sep = XPLMGetDirectorySeparator();
    if ( inifiles[(long)current_inifile] == "smoke.ini" ) {
      std::string tmp;
      int r = ac_directory(tmp);
      if ( r == 1 ) {
	G.read_prefs( tmp );
      } else {
	G.read_prefs( the_dir + sep + inifiles[(long)current_inifile] );
      }
    }
    G.read_smoke_prefs( the_dir + sep + "ai.ini", G.ai_smoke );
    G.read_smoke_prefs( the_dir + sep + "wt.ini", G.wt_smoke );
    lg.xplm( "wt.ini rgb "+std::to_string(G.wt_smoke.rgb_right[0])+"\n" );
  }
  
  if ( inMessage == XPLM_MSG_LIVERY_LOADED ) { // 108
    //rescan = true; // this disturbs the loading above in 102?
  }
  
  if ( inMessage == XPLM_MSG_AIRPORT_LOADED ) { // 103
    smokemode = false;
    aismokemode = false;
    wtsmokemode = false;
    rescan = true;
    XPLMCheckMenuItem(myMenu, MENU_TOGGLE, xplm_Menu_Unchecked);
    XPLMCheckMenuItem(myMenu, MENU_AI_TOGGLE, xplm_Menu_Unchecked);
    XPLMCheckMenuItem(myMenu, MENU_WT_TOGGLE, xplm_Menu_Unchecked);
  }
  
}

/**********************************************************************/

int SmokeCommandHandler(XPLMCommandRef inCommand, XPLMCommandPhase inPhase, void *inRefcon) {
  if (inPhase == 0) {
    smokemode = ! smokemode;
    // copied from menu handler:
    if ( smokemode ) {
      XPLMCheckMenuItem(myMenu, MENU_TOGGLE, xplm_Menu_Checked);
    } else {
      XPLMCheckMenuItem(myMenu, MENU_TOGGLE, xplm_Menu_Unchecked);
    }
  }
  return 0;
}

int AISmokeCommandHandler(XPLMCommandRef inCommand, XPLMCommandPhase inPhase, void *inRefcon) {
  if (inPhase == 0) {
    aismokemode = ! aismokemode;
    // copied from menu handler:
    if ( aismokemode ) {
      XPLMCheckMenuItem(myMenu, MENU_AI_TOGGLE, xplm_Menu_Checked);
    } else {
      XPLMCheckMenuItem(myMenu, MENU_AI_TOGGLE, xplm_Menu_Unchecked);
    }
  }
  return 0;
}

int WTSmokeCommandHandler(XPLMCommandRef inCommand, XPLMCommandPhase inPhase, void *inRefcon) {
  if (inPhase == 0) {
    wtsmokemode = ! wtsmokemode;
    // copied from menu handler:
    if ( wtsmokemode ) {
      XPLMCheckMenuItem(myMenu, MENU_WT_TOGGLE, xplm_Menu_Checked);
    } else {
      XPLMCheckMenuItem(myMenu, MENU_WT_TOGGLE, xplm_Menu_Unchecked);
    }
  }
  return 0;
}

int EulerCommandHandler(XPLMCommandRef inCommand, XPLMCommandPhase inPhase, void *inRefcon) {
  if (inPhase == 0) {
    if ( headmode == -1 ) { // we enable, look for closes, or first, or ...
      headmode = closest_mp();
      findapmode = -1;
    } else {
      headmode = -1; // disable head tracking
      dr_ph_psi = 0; // restore
      dr_ph_the = 0;
      dr_ph_phi = 0; 
    }
  }
  return 0;
}

int FindAPCommandHandler(XPLMCommandRef inCommand, XPLMCommandPhase inPhase, void *inRefcon) {
  if (inPhase == 0) {
    if ( findapmode == -1 ) { // we enable, look for closes, or first, or ...
      float lat, lon;
      G.get_nearest_ap( dr_plane_lat, dr_plane_lon, lat, lon ); // also in G.nearest_ap_lat/lon
      findapmode = 1;
      headmode = -1;
    } else {
      findapmode = -1; // disable head tracking
      dr_ph_psi = 0; // restore
      dr_ph_the = 0;
      dr_ph_phi = 0; 
    }
  }
  return 0;
}

void MyMenuHandlerCallback( void *inMenuRef, void *inItemRef) {
  (void)inMenuRef;

  if ( (intptr_t)inItemRef == MENU_TOGGLE ) {
    smokemode = ! smokemode;
    if ( smokemode ) {
      XPLMCheckMenuItem(myMenu, MENU_TOGGLE, xplm_Menu_Checked);
    } else {
      XPLMCheckMenuItem(myMenu, MENU_TOGGLE, xplm_Menu_Unchecked);
    }	
  }
  if ( (intptr_t)inItemRef == MENU_AI_TOGGLE ) {
    aismokemode = ! aismokemode;
    if ( aismokemode ) {
      XPLMCheckMenuItem(myMenu, MENU_AI_TOGGLE, xplm_Menu_Checked);
    } else {
      XPLMCheckMenuItem(myMenu, MENU_AI_TOGGLE, xplm_Menu_Unchecked);
    }	
  }
  if ( (intptr_t)inItemRef == MENU_WT_TOGGLE ) {
    wtsmokemode = ! wtsmokemode;
    if ( wtsmokemode ) {
      XPLMCheckMenuItem(myMenu, MENU_WT_TOGGLE, xplm_Menu_Checked);
    } else {
      XPLMCheckMenuItem(myMenu, MENU_WT_TOGGLE, xplm_Menu_Unchecked);
    }	
  }
  if ( (intptr_t)inItemRef == MENU_RELOAD ) {
    /*
    lg.xplm("Reloading config file.\n");
    delete_textures();
    load_textures(); // doesn't work, needs unload first?
    G.read_prefs();
    rescan = true;
    lg.xplm("Reloaded config file.\n");
    */
    // The above should be rescan-directory, and then a list with ini files
    //
    // Remove old entries
    XPLMCheckMenuItem(myMenu, (long)(current_inifile+5), xplm_Menu_NoCheck); // OFFSET 5
    for( int i = 0; i < inifiles.size(); i++ ) {
      XPLMRemoveMenuItem(myMenu, 5); // always 5... because they move after one delete...
    }
    //    
    inifiles.clear();
    m_particleList.clear();
    current_inifile = 0;
    scan_directory(the_dir, inifiles);
    for( int i = 0; i < inifiles.size(); i++ ) {
      XPLMAppendMenuItem(myMenu, inifiles[i].c_str(),  (void*)((uintptr_t)(i+5)), 1); // OFFSET 5
      if ( inifiles[i] == "smoke.ini" ) {
	current_inifile = i;
      }
      // PJB TODO look for a wt.ini and ai.ini and save those
    }
    XPLMCheckMenuItem(myMenu, ((long)current_inifile+5), xplm_Menu_Checked); // OFFSET 5
        
    const char *sep = XPLMGetDirectorySeparator();
    if ( inifiles[(long)current_inifile] != "smoke.ini" ) {
      G.read_prefs( the_dir + sep + inifiles[(long)current_inifile] );
    } else {
      std::string tmp;
      int r = ac_directory(tmp);
      if ( r == 1 ) {
	G.read_prefs( tmp );
      } else {
	G.read_prefs( the_dir + sep + inifiles[(long)current_inifile] );
      }
    }
    G.read_smoke_prefs( the_dir + sep + "ai.ini", G.ai_smoke );
    G.read_smoke_prefs( the_dir + sep + "wt.ini", G.wt_smoke );
    lg.xplm( "wt.ini rgb "+std::to_string(G.wt_smoke.rgb_right[0])+"\n" );
  } // RELOAD MENU
  
  // from 5 onwards are config files, load them here on the fly, like G.read_prefs, with
  // inifiles[((intptr_t)inItemRef-3] as argument
  //
  lg.xplm("MENU="+std::to_string((intptr_t)inItemRef)+"/"+std::to_string(inifiles.size())+"\n");
  if ( (intptr_t)inItemRef > 4 ) { // OFFSET
    XPLMCheckMenuItem(myMenu, (long)(current_inifile+5), xplm_Menu_NoCheck); // uncheck curent one, OFFSET 4
    current_inifile = (intptr_t)inItemRef-5; // OFFSET
    lg.xplm("CURRENT="+std::to_string(current_inifile)+"/"+std::to_string(inifiles.size())+"\n");
    const char *sep = XPLMGetDirectorySeparator();
    if ( inifiles[(long)current_inifile] != "smoke.ini" ) {
      G.read_prefs( the_dir + sep + inifiles[(long)current_inifile] );
    } else {
      std::string tmp;
      int r = ac_directory(tmp);
      if ( r == 1 ) {
	G.read_prefs( tmp );
      } else {
	G.read_prefs( the_dir + sep + inifiles[(long)current_inifile] );
      }
    }
    G.read_smoke_prefs( the_dir + sep + "ai.ini", G.ai_smoke );
    G.read_smoke_prefs( the_dir + sep + "wt.ini", G.wt_smoke );
    lg.xplm( "wt.ini rgb "+std::to_string(G.wt_smoke.rgb_right[0])+"\n" );

    XPLMCheckMenuItem(myMenu, (intptr_t)inItemRef, xplm_Menu_Checked);

    rescan = true;
  }
}

// Do the real initialisation when ready then change to the "real" flight loop
// (not really used here anymore)
float DeferredInitNewAircraftFLCB(float elapsedMe, float elapsedSim, int counter, void * refcon) {
  static int MyProgramFLCBStartUpFlag = 0;
  //lg.xplm("DeferredInitNewAircraftFLCB finished\n");
  XPLMRegisterFlightLoopCallback(MyFlightLoopCallback, G.smoke_refresh, NULL);
  XPLMRegisterFlightLoopCallback(MyHeadLoopCallback, 1.0, NULL);
  XPLMRegisterFlightLoopCallback(MyFindAPLoopCallback, 1.0, NULL);

  lg.xplm("G.init_mp_datarefs()\n");
  G.init_mp_datarefs();
  
  return 0; // Returning 0 stops DeferredInitFLCB from being looped again.
}

float MyFlightLoopCallback(float inElapsedSinceLastCall, float inElapsedTimeSinceLastFlightLoop, int inCounter, void *inRefcon) {

  if ( inElapsedSinceLastCall > 1.0 ) { //we've been in the background or something
    return G.smoke_refresh;
  }
  if ( dr_sim_tfts == last_tfts ) { // this probably stops replay-moke, but allows in map?
    //return G.smoke_refresh;
  }
  last_tfts = dr_sim_tfts;
	
  if ( dr_sim_paused == 1 ) { // or if in map view?
    return G.smoke_refresh;
  }
  
  int    sta     = dr_sim_speed;
  double elapsed = inElapsedSinceLastCall * sta;

  bool emit = false;
  if ( m_particleList.size() < G.max_bols ) {
    emit = true;
  } else {
    //lg.xplm( "No particles available.\n" );
  }

  if ( rescan ) {
    scan_engines();
    rescan = false;
  }
  float smoke_emit_rate = G.plane_smoke_emit_rate;
  float throttle_factor = 0.25 + (dr_throttle_ratio_all*0.75);
  if ( G.plane_smoke_throttle > 0 ) {
    //smoke_emit_rate = G.plane_smoke_emit_rate * (0.5 + (dr_throttle_ratio_all/2.0));
    smoke_emit_rate = G.plane_smoke_emit_rate * throttle_factor;
    //lg.xplm("smoke_emit_rate="+rounded(smoke_emit_rate)+"\n");
  }
  if ( (smoke_emit_rate > 0) && smokemode ) { // arbitrary minimum speed?
    if ( m_particleList.size() < G.max_bols ) {
      float smoke_flow = G.plane_smoke_flow;
      if ( G.plane_smoke_throttle > 0 ) {
	//smoke_flow = G.plane_smoke_flow * (0.5 + (dr_throttle_ratio_all/2.0));
	smoke_flow = G.plane_smoke_flow * throttle_factor;
	//lg.xplm("smoke_flow="+rounded(smoke_flow)+"\n");
      }
      if ( G.p_rand() <= smoke_flow ) {    // "thickness", "puffiness"
	// precalc
	float phi = dr_plane_phi*(M_PI/180);
	float psi = dr_plane_psi*(M_PI/180);
	float the = dr_plane_the*(M_PI/180);
	float alt = 10;// dr_plane_y_agl;
	float d_agl = dr_plane_y_agl - dr_plane_ly; // delta agl to calc ground level, so we don't need probe in Bol.h
	//float alt = dr_plane_elevation; //is a double, but precision not so important
	float x_wrl;
	float y_wrl;
	float z_wrl;
	//
	// calculate from plane_x_prev &c...
	float x_wrl_prev;
	float y_wrl_prev;
	float z_wrl_prev;
	//
	float wind_vel;
	float wind_dir;
	G.interpolate_wind(alt, wind_vel, wind_dir);
	//lg.xplm("wind0 "+std::to_string(wind_dir)+","+std::to_string(wind_vel)+"\n");
	//
	float rd = G.random_range(-12,12); // random direction change, make parameter?
	float wind_s =  sin( (wind_dir+(180.0+rd)) * (float)(M_PI/180.0) );
	float wind_c = -cos( (wind_dir-(180.0-rd)) * (float)(M_PI/180.0) );
	//
	//lg.xplm("wind0 "+std::to_string(wind_s)+","+std::to_string(wind_c)+"\n");
	//lg.xplm("wind1 "+std::to_string(dr_wind_now_x_msc)+","+std::to_string(dr_wind_now_z_msc)+"\n");
	//
	// dr_wind_direction_degt
	// dr_wind_speed_kt
	wind_dir = dr_wind_direction_degt;
	wind_vel = dr_wind_speed_kt * 0.5144; // knots to m/s
	// a wind_vel of 0 looks strange... if no wind, make wind
	if ( wind_vel < 1.0 ) {
	  wind_vel = 1.0;
	}
	/*
	wind_s = dr_wind_now_x_msc + pn_rand()*0.10; 
	wind_c = dr_wind_now_z_msc + pn_rand()*0.10;
	*/
	// height
	// plane is at y_agl (height above ground), and we have y-coordinate.
	// zero height is at this point: y-coordinate - y_agl
	//
	/*
	  pp.psi   = dr_psi;
	  pp.pitch = dr_theta;
	  pp.roll  = dr_phi;
	  conversion(x_plane, y_plane, z_plane, // = source location in airplane coordinates.  
	  phi, psi, the, 
	  local_x, local_y, local_z,  //plane's location in the world 
	  x_wrl, y_wrl, z_wrl) ;
	*/
	// plane_smoke_emitters = 0, 5.5, 0.25
	// plane_smoke_emitters = 4, 5, 2
	// the emitters
	////int emits = G.plane_smoke_emit_rate * elapsed * G.random_variation(30);
	int emits = smoke_emit_rate * elapsed * G.random_variation(30);
	//lg.xplm( "A emits="+rounded(emits)+" emits_subcounter="+rounded(emit_subcounter)+" elap="+rounded(elapsed)+"\n" );
	//  maybe do this always? even for larger numbers?
	if ( emits < 1 ) {
	  //lg.xplm( "emits_subcounter="+rounded(emit_subcounter)+"\n" );
	  emit_subcounter += (float)smoke_emit_rate * elapsed;
	  if ( emit_subcounter >= 1 ) { //G.plane_smoke_emit_rate ) { // was 1.0 (whould be 1/G.plane_smoke_emit_rate ?)
	    emits = 1;
	    emit_subcounter = 0.0;
	  }
	}
	//lg.xplm( "B emits="+rounded(emits)+" emits_subcounter="+rounded(emit_subcounter)+"\n" );
	//
	if ( emits > 0 ) {
	  float smoke_emit_velocity = G.plane_smoke_emit_velocity;
	  if ( G.plane_smoke_throttle > 0 ) {
	    smoke_emit_velocity = G.plane_smoke_emit_velocity * throttle_factor;
	  }

	for (int i = 0; i < G.emit_num; i++ ) {
	  // rate should really be per sec (e.g. 100, then multiply with elapsed time).
	  // 30 fps = 1/30 = 0.0333, to get 8 puffs, we need a setting of 240.
	  // but the distance travelled should then also be taken into account...
	  // we have sim/flightmodel/position/true_airspeed available, a distance
	  // could be added to the z offset?
	  // easiest to save last x,y,z pos, then interpolate between? x1 + (x0-x1)*p_rand(),  etc
	  
	  // Calculate position from plane coordinates to world coordinates.
	  conversion( G.emit_x[i], G.emit_y[i], G.emit_z[i], //+lvz, 
		      phi, psi, the,
		      dr_plane_lx, dr_plane_ly, dr_plane_lz, // for +lvz we need this static outside loop?
		      x_wrl, y_wrl, z_wrl);
	  // Previous position
	  conversion( G.emit_x[i], G.emit_y[i], G.emit_z[i], //+lvz, 
		      phi, psi, the,
		      plane_x_prev, plane_y_prev, plane_z_prev,
		      x_wrl_prev, y_wrl_prev, z_wrl_prev);
	  float x_wrl_inc = x_wrl_prev - x_wrl; // The position delta
	  float y_wrl_inc = y_wrl_prev - y_wrl;
	  float z_wrl_inc = z_wrl_prev - z_wrl;

	  //float xyz_dist = sqrt((x_wrl_inc*x_wrl_inc)+(y_wrl_inc*y_wrl_inc)+(z_wrl_inc*z_wrl_inc));
	  //lg.xplm("xyz_dist="+rounded(xyz_dist)+"\n" );
	  
	  int speed_factor = 1;
	  /*
	  if ( (fabs(x_wrl_inc) > 0.2) || (fabs(z_wrl_inc) > 0.2) ) {
	    speed_factor = 4;
	  }
	  if ( (fabs(x_wrl_inc) > 0.5) || (fabs(z_wrl_inc) > 0.5) ) {
	    speed_factor = 8;
	  }
	  if ( (fabs(x_wrl_inc) > 0.8) || (fabs(z_wrl_inc) > 0.8) ) {
	    speed_factor = 12;
	  }
	  */
	  
	  //lg.xplm("x_wrl="+rounded(x_wrl)+", y_wrl="+rounded(y_wrl)+", z_wrl="+rounded(z_wrl)+"\n" );
	  //lg.xplm("x_prv="+rounded(x_wrl_prev)+", y_prv="+rounded(y_wrl_prev)+", z_prv="+rounded(z_wrl_prev)+"\n" );
	  //lg.xplm("x_inc="+rounded(x_wrl_inc)+", y_inc="+rounded(y_wrl_inc)+", z_inc="+rounded(z_wrl_inc)+"\n" );
	  
	  //for( int x = 0; x < G.plane_smoke_emit_rate; x += 1) { // "emit_rate"
	  for( int x = 0; x < emits*speed_factor; x += 1) { // "emit_rate"
	    SmokeBol *b = new SmokeBol();
	    /*
	    // Calculate position from plane coordinates to world coordinates.
	    conversion( G.emit_x[i], G.emit_y[i], G.emit_z[i], //+lvz, 
			phi, psi, the,
			dr_plane_lx, dr_plane_ly, dr_plane_lz, // for +lvz we need this static outside loop?
			x_wrl, y_wrl, z_wrl);
	    // Previous position
	    conversion( G.emit_x[i], G.emit_y[i], G.emit_z[i], //+lvz, 
			phi, psi, the,
			plane_x_prev, plane_y_prev, plane_z_prev,
			x_wrl_prev, y_wrl_prev, z_wrl_prev);
	    */
	    // Position based on movement
	    // PJB: should the "+ G.pn_rand() * G.plane_smoke_movement" be here? Without it, it starts more
	    // compact...
	    //
	    // If we have previous x_wrl, y_wrl and z_wrl, we pick a position between
	    // the new x_wrl, and the previous ones.
#ifdef STARTMOVEMENT
	    b->x = x_wrl + G.pn_rand() * G.plane_smoke_movement; //G.pn_rand()/2.0; //pp.lx + (s * dist_offset) + G.pn_rand();
	    b->y = y_wrl; //pp.ly; // altitude
	    b->z = z_wrl + G.pn_rand() * G.plane_smoke_movement;//G.pn_rand()/2.0; //pp.lz + (c * dist_offset) + G.pn_rand();
#else
	    // Interpolate between the current and last position
	    b->x = x_wrl;
	    b->y = y_wrl;
	    b->z = z_wrl;
	    //
	    //
	    /*
	    lg.xplm("i="+rounded(i)+", x="+rounded(x)+"\n");
	    lg.xplm("x_wrl="+rounded(x_wrl)+", y_wrl="+rounded(y_wrl)+", z_wrl="+rounded(z_wrl)+"\n" );
	    lg.xplm("x_prv="+rounded(x_wrl_prev)+", y_prv="+rounded(y_wrl_prev)+", z_prv="+rounded(z_wrl_prev)+"\n" );
	    lg.xplm("x_inc="+rounded(x_wrl_inc)+", y_inc="+rounded(y_wrl_inc)+", z_inc="+rounded(z_wrl_inc)+"\n" );
	    */
	    //
	    // Now we seem spread out more at higher speeds...why?
	    float foo = G.p_rand();
	    b->x = x_wrl + x_wrl_inc * foo;//G.p_rand();
	    b->y = y_wrl + y_wrl_inc * foo;//G.p_rand();
	    b->z = z_wrl + z_wrl_inc * foo;//G.p_rand();
#endif
	    // give the particle a random velocity (plus planes with decelleration?)
	    //
	    // -> or always 180 degs from heading? up/down?
	    float emit_dir = dr_plane_psi + 180;
	    float emit_dir_s =  sin( emit_dir * (float)(M_PI/180.0) );
	    float emit_dir_c = -cos( emit_dir * (float)(M_PI/180.0) );
	    // emit velocity... ?
	    //float emit_vel = G.plane_smoke_emit_velocity * p_rand()+0.5; // or random_variation()?
	    float emit_vel = smoke_emit_velocity * G.random_variation(30); // or random_variation()?
	    // emit_val throttle related? 
	    // G.random_variation( 20,30 ) => -20 and plus 30% ?
	    //                "* G.plane_smoke_movement" is needed to keep e.g black smoke together
	    b->vx = emit_dir_s * G.plane_smoke_movement * emit_vel; // with high vel maybe not "* emit_vel"?
	    b->vz = emit_dir_c * G.plane_smoke_movement * emit_vel; // ^^no, it seems needed
	    // Add plane speed? (subtract, other direction)
	    //b->vx -= dr_plane_lvx; (doesn't seem to work properly)
	    //b->vz -= dr_plane_lvz;
	    //
	    //b->vx = G.pn_rand() * G.plane_smoke_movement;   
	    b->vy = G.plane_smoke_start_vy * G.random_variation(30); //G.p_rand()  * G.plane_smoke_movement;
	    //b->vz = G.pn_rand() * G.plane_smoke_movement;
	    //
	    float min_size = G.plane_smoke_size_from * G.random_variation(30);
	    float max_size = G.plane_smoke_size_to   * G.random_variation(30);
	    float growth   = G.plane_smoke_size_grow * G.random_variation(30);
	    //    red = 1.0*p_rand() ?
	    float fudge = pn_rand() * G.plane_smoke_colour_var;
	    float re, gr, bl;
	    if ( G.emit_x[i] > 0 ) {
	      re = G.plane_smoke_rgb_right[0] + fudge;
	      gr = G.plane_smoke_rgb_right[1] + fudge;
	      bl = G.plane_smoke_rgb_right[2] + fudge;
	    } else {
	      re = G.plane_smoke_rgb_left[0] + fudge;
	      gr = G.plane_smoke_rgb_left[1] + fudge;
	      bl = G.plane_smoke_rgb_left[2] + fudge;
	    }
	    //
	    float dd = G.plane_smoke_draw_delay;
	    float tf = G.plane_smoke_transparency_from;
	    float tt = G.plane_smoke_transparency_to;
	    //
	    b->set_linger( G.plane_smoke_linger_min + (G.plane_smoke_linger_time*p_rand()), G.plane_smoke_linger_grow);
	    b->set_wind( wind_s, wind_c, wind_vel*G.plane_smoke_wind_factor );
	    //b->set_wind( wind_s*G.plane_smoke_wind_factor, wind_c*G.plane_smoke_wind_factor );
	    b->set_sizes(min_size, max_size, growth,  re, gr, bl,   dd,    tf, tt);
	    b->set_target_vy(G.plane_smoke_target_vy);
	    //
	    // y_offset interfers? random spawn pos interfers? If vy is downwards...
	    b->ground_zero = dr_plane_ly - dr_plane_y_agl - 0.20; // - (max_size/2.0); // stick up a little
	    //lg.xplm( "dr_plane_ly="+rounded(dr_plane_ly)+", dr_plane_y_agl="+rounded(dr_plane_y_agl)+"\n" );
	    //
	    m_particleList.push_back(b);
	  } // x
	} // i
	} // emits > 0
      }
    }
  }
  // we probably have moved already...
  plane_x_prev = dr_plane_lx;
  plane_y_prev = dr_plane_ly;
  plane_z_prev = dr_plane_lz;

  /*
  if ( G.p_rand() > 0.9 ) {
    RisingSmokeBol *b = new RisingSmokeBol();
    float c = G.p_rand();
    //           sizes      growth  colours   delay    transparency
    b->set_sizes(1.0, 30.0, 0.5,    c, c, c,   0.1,    0.8, 0.0);
    b->x = dr_plane_lx;
    b->y = dr_plane_ly;
    //b->vy = 1;
    b->z = dr_plane_lz;
    m_particleList.push_back(b);
  }
  */

  if ( findapmode > 0 ) { // when lookat mode is enabled
    if ( G.p_rand() > 0.9 ) {
      RisingSmokeBol *b = new RisingSmokeBol();
      float c = G.p_rand();
      //           sizes      growth  colours   delay    transparency
      b->set_sizes(1.0, 30.0, 0.5,    c, c, c,   0.1,    0.8, 0.0);
      b->x = G.nearest_ap_x;
      b->y = G.nearest_ap_y + 100; // 100 m above airport (otherise ON rwy)
      //b->vy = 1;
      b->z = G.nearest_ap_z;
      m_particleList.push_back(b);
    }
  }
  
  // Put a smoke puff at AI plane positions.
  if ( aismokemode ) {
    int num_mp = G.get_mp_positions();
    for (int i = 0; i < num_mp; i++ ) {
      //std::string m =" POS AI"+std::to_string(i)+" :"+std::to_string(G.mp_xs[i])+", "+std::to_string(G.mp_ys[i])+", "+std::to_string(G.mp_zs[i])+std::to_string(G.mp_vxs[i])+", "+std::to_string(G.mp_vys[i])+", "+std::to_string(G.mp_vzs[i])+"\n";
      //lg.xplm(m);
      if ( G.p_rand() < G.ai_smoke.flow ) {
	if ( (fabs(G.mp_vxs[i]) > 1) || (fabs(G.mp_vzs[i]) > 1) ) { 
	  SmokeBol *b = new SmokeBol();
	  float fudge = pn_rand() * G.ai_smoke.colour_var;
	  float re, gr, bl;
	  re = G.ai_smoke.rgb_right[0] + fudge;
	  gr = G.ai_smoke.rgb_right[1] + fudge;
	  bl = G.ai_smoke.rgb_right[2] + fudge;
	  b->set_sizes( G.ai_smoke.size_from, G.ai_smoke.size_to, G.ai_smoke.size_grow,
			re, gr, bl,
			G.ai_smoke.draw_delay,
			G.ai_smoke.transparency_from, G.ai_smoke.transparency_to );
	  b->set_linger( G.ai_smoke.linger_min + (G.ai_smoke.linger_time*p_rand()), G.ai_smoke.linger_grow);
	  b->x = G.mp_xs[i] + G.pn_rand();
	  b->y = G.mp_ys[i] + G.pn_rand();
	  b->z = G.mp_zs[i] + G.pn_rand();
	  //
	  b->vx = G.pn_rand() * 0.1;
	  b->vy = G.pn_rand() * 0.1;
	  b->vz = G.pn_rand() * 0.1;
	  //
	  m_particleList.push_back(b);
	  //std::string m =" POS AI"+std::to_string(i)+" :"+std::to_string(G.mp_xs[i])+", "+std::to_string(G.mp_ys[i])+", "+std::to_string(G.mp_zs[i])+"\n";
	  //lg.xplm(m);
	}
      }
    }
  }

  if ( wtsmokemode ) {
    int foo = G.get_cjs();
    for( int i = 0; i < foo; i++ ) {
      // dr_plane_lat/lon
      //float lat = static_cast<float>(plane_lat);
      double l_x, l_y, l_z;
      float lat = G.cjs_ac_lat[i];
      float lon = G.cjs_ac_lon[i];
      float asl = G.cjs_ac_asl[i];
      float spd = G.cjs_speed_kias[i];
      //float dt = distanceto(dr_plane_lat, dr_plane_lon, lat, lon); // Not used
      XPLMWorldToLocal( lat, lon, asl/3.28084, &l_x, &l_y, &l_z );
      if ( spd > 0 ) {
	if ( spd < 10000.0 ) { // some filtering to remove rubbish. PJB FIXME
	  if ( G.p_rand() < G.wt_smoke.flow ) {
	    SmokeBol *b = new SmokeBol();
	    float fudge = 0;//pn_rand() * G.wt_smoke.colour_var; // PJB fudge to big? FIXME
	    float re, gr, bl;
	    re = G.wt_smoke.rgb_right[0] + fudge;
	    gr = G.wt_smoke.rgb_right[1] + fudge;
	    bl = G.wt_smoke.rgb_right[2] + fudge;
	    //std::string m =" RGB: "+std::to_string(re)+", "+std::to_string(gr)+", "+std::to_string(bl)+"\n";
	    //lg.xplm(m);
	    b->set_sizes( G.wt_smoke.size_from, G.wt_smoke.size_to, G.wt_smoke.size_grow,
			  re, gr, bl,
			  G.wt_smoke.draw_delay,
			  G.wt_smoke.transparency_from, G.wt_smoke.transparency_to );
	    b->set_linger( G.wt_smoke.linger_min + (G.wt_smoke.linger_time*p_rand()), G.wt_smoke.linger_grow);
	    //b->set_wind( wind_s, wind_c, wind_vel*G.plane_smoke_wind_factor );
	    b->x = l_x + G.pn_rand();
	    b->y = l_y + G.pn_rand();
	    b->z = l_z + G.pn_rand();
	    b->vx = G.pn_rand() * 0.1;
	    b->vy = G.pn_rand() * 0.1;
	    b->vz = G.pn_rand() * 0.1;
	    m_particleList.push_back(b);
	  }
	  //
	  //lg.xplm("CJS: "+std::to_string(i)+", "+std::to_string(lat)+", "+std::to_string(lon)+", "+std::to_string(asl)+", "+std::to_string(spd)+", "+std::to_string(dt)+"\n");
	}
      }
    }
  }

  /*
  float fps = 1.0 / inElapsedSinceLastCall;
  lg.xplm("PARTICLES: "+std::to_string((int)m_particleList.size())+" FPS="+rounded(fps)+"\n");
  */
  for (xcb_i = m_particleList.begin(); xcb_i != m_particleList.end();) {
    CBol *b = *xcb_i;
    int res = b->update(elapsed); // if status 0, remove?
    if ( res == 0 ) {
      xcb_i = m_particleList.erase(xcb_i);
    } else {
      ++xcb_i;
    }
  }

  // Return time interval till next call
  return G.smoke_refresh;
}

// Separate flightloop for head tracking 
float MyHeadLoopCallback(float inElapsedSinceLastCall, float inElapsedTimeSinceLastFlightLoop, int inCounter, void *inRefcon) {
  if ( headmode > -1 ) {
    euler( headmode );
  }
  return G.smoke_refresh; //1.0;
}

// Seperate flightloop for head tracking 
float MyFindAPLoopCallback(float inElapsedSinceLastCall, float inElapsedTimeSinceLastFlightLoop, int inCounter, void *inRefcon) {
  if (findapmode > -1 ) {
    move_head_towards(G.nearest_ap_lat, G.nearest_ap_lon, 100); // airport
  }
  return G.smoke_refresh; //1.0;
}


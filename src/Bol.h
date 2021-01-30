
#ifndef _BOL_H
#define _BOL_H

#include <string>
#include <vector>
#include <random>
#include "lodepng.h"
#include "Log.h"

#define PI 3.14159265358979323846

// ----------------------------------------------------------------------------
// Class
// ----------------------------------------------------------------------------

/*
sim/graphics/view/
plane_render_type	int	1030+	no	enum	What kind of rendering is being done during the xplm_Phase_Airplanes callback. solid = 1, blended/alpha = 2 (0 = N/A, outside draw callback)

 */
namespace SMOKE {
 
  float p_rand() { // positive random 0..1 copy of Global version....
    return float(std::rand()) / float(RAND_MAX);
  }
  float pn_rand() { // positive-negative random -1..1
    return (float(std::rand())-float(std::rand())) / float(RAND_MAX);
  }

  static XPLMDataRef plane_render_type = XPLMFindDataRef("sim/graphics/view/plane_render_type"); // int 10.30+
  static XPLMDataRef world_render_type = XPLMFindDataRef("sim/graphics/view/world_render_type"); // int 10.30+
  static XPLMDataRef hdr_on = XPLMFindDataRef("sim/graphics/settings/HDR_on"); // int, 10.25+
  //DataRef<float>  dr_sun_pitch("sim/graphics/scenery/sun_pitch_degrees");
  static XPLMDataRef sun_pitch = XPLMFindDataRef("sim/graphics/scenery/sun_pitch_degrees"); // int, 10.25+
  static XPLMDataRef gModelviewMatrixRef = XPLMFindDataRef("sim/graphics/view/modelview_matrix");
  
  // -- own local version
  static XPLMProbeRef hProbeBol = XPLMCreateProbe(xplm_ProbeY); // the includes?
  float height(double x, double y, double z) {
    XPLMProbeInfo_t info = { 0 };
    info.structSize = sizeof(info);  
    // If we have a hit then return Y coordinate
    if (XPLMProbeTerrainXYZ( hProbeBol, x, y, z, &info) == xplm_ProbeHitTerrain) {
      return info.locationY;
    }
    return -1.0;
  }
  // --


  void camera_directions(float * out_rgt,// Any can be NULL
				float * out_up ,
				float * out_look) {
    float m[16];
    glGetFloatv(GL_MODELVIEW_MATRIX, m);
    
    // Roughly speaking, a modelview matrix is made up more or less like this:
    // [ EyeX_x EyeX_y EyeX_z    a
    //   EyeY_x EyeY_y EyeY_z    b
    //   EyeZ_x EyeZ_y EyeZ_z    c
    //   um don't look down here ]
    // where a, b, c are translations in _eye_ space.  (For this reason, a,b,c is not
    // the camera location - sorry!)
    
    if (out_rgt) {
      out_rgt[0] = m[0];
      out_rgt[1] = m[4];
      out_rgt[2] = m[8];
    }
    if (out_up) {
      out_up[0] = m[1];
      out_up[1] = m[5];
      out_up[2] = m[9];
    }
    if (out_look) {
      out_up[0] = m[2];
      out_up[1] = m[6];
      out_up[2] = m[10];
    }
  }
  
  XPLMTextureID circle_texture[1];
  std::vector<unsigned char> image;
  
  void load_textures() {
    std::string filename = "Resources/plugins/xsmoke/smoke.png";
    
    unsigned w, h;
    unsigned error = lodepng::decode(image, w, h, filename.c_str());
    if ( error ) {
      XPLMDebugString(lodepng_error_text(error));
    } else {
      lg.xplm( "Loaded texture "+filename+"\n");
      std::string x = "w="+std::to_string(w)+",h="+std::to_string(h)+"\n";
      lg.xplm( x );
    }
    XPLMGenerateTextureNumbers( &circle_texture[0], 1 ); 
    XPLMBindTexture2d( circle_texture[0], 0 );
    glTexImage2D( GL_TEXTURE_2D, 0, GL_RGBA, w, h, 0, GL_RGBA, GL_UNSIGNED_BYTE, &image[0] ); // image!

    glTexParameteri( GL_TEXTURE_2D,GL_TEXTURE_MIN_FILTER,GL_LINEAR );
    glTexParameteri( GL_TEXTURE_2D,GL_TEXTURE_MAG_FILTER,GL_LINEAR );
  }

  void delete_textures() {
    glDeleteTextures(1, (const GLuint *)&circle_texture[0]);
  }
  
static int bol_counter = 0;
  
  class CBol {
    
  public:
    std::string name;
    std::string num;
    size_t idx;

    int status;
    float x, y, z, r;
    float vx = 0.0;
    float vy = 0.0;
    float vz = 0.0;
    // accel, time, etc
    float age;
    float tp = 1.0; // initial transparency

    int typ = 0;  
    float min_size = 0.01;
    float max_size = 0.20;
    float growth   = 0.05; //growth in m/s
    //
    float max_age = (max_size-min_size)/growth; // 0.19 / 0.05 = 3.8
    float trans_from = 0.3;
    float trans_to = 0.1;
    //float trans_growth = growth/(max_size-min_size)*(trans_from-trans_to); // 0.05 / (0.19)*(0.2) = 0.05 / 0.038 = 1.32
    float trans_growth = (trans_from-trans_to) / max_age;
    //
    float draw_delay = 0.0;

    float linger = 10.0;
    float linger_grow = 1.0;
    
    // Accel
    float vx_a = 0.0;
    float vy_a = 0.0;
    float vz_a = 0.0;

    // colour?
    float smoke_rgb[3];

    float wind_x = 1.0;
    float wind_z = 1.0;

    // zero height, arbitrary really low
    float ground_zero = -30000.0;

    float mModelView[16], mProjection[16], m[16];
    
    // Constructor.
    CBol() {
      status = 1;
      age = 0.0;
      smoke_rgb[0] = 0.8;
      smoke_rgb[1] = 0.8;
      smoke_rgb[2] = 0.8;
    };
    // Destructor.
    ~CBol() {} ;
    
    virtual int update(float t) {
      return 0;
    }

  };

class Bol : public CBol {
 public:
  int typ = 0;

  int update(float delta_t) {
    if ( status == 0 ) {
      return 0;
    }
    if ( status == 1 ) {
      r += 10.0 * delta_t;
      if ( r >= 5.0 ) {
	status = 2;
      }
    } else if ( status == 2 ) {
      r -= 5.0 * delta_t;
      if ( r <= 0.1 ) {
	status = 0; // die
      }
    }
    tp -= 0.01; //need a nice function for this
    age += delta_t;
    return status;
  }
};

/*
  Extra parameter "linger" in grow, which is extra the max age at last setting (don't make transp. 0 in that case)
  smoke_size = 1.0, 3.0, 0.1, 10.0
                              ____
 */

class SmokeBol : public CBol {
 public:
  // Constructor.
  SmokeBol() {
    status = 1;
    age = 0.0;
    tp = trans_from;
    r = min_size; // start size, from 0.1 to 5 is 4.9, 1m/s is 4.9 seconds maxage
    linger = 60.0 * p_rand();
  };

  void set_wind(float wx, float wz, float wv) { // wind_vel
    wind_x = wx*(wv+0.01)*(p_rand()+0.5);
    wind_z = wz*(wv+0.01)*(p_rand()+0.5);
    // recalc accel? need age and max_age
    vx_a = (wind_x - vx) / max_age; // in max_age, then stay at speed in linger!
    vz_a = (wind_z - vz) / max_age; // in max_age, then stay at speed in linger!
  }

  // if we take the x-plane supplied wind vectors
  void set_wind(float wx, float wz) { 
    wind_x = wx;
    wind_z = wz;
    // recalc accel? need age and max_age
    vx_a = (wind_x - vx) / max_age; // in max_age, then stay at speed in linger!
    vz_a = (wind_z - vz) / max_age; // in max_age, then stay at speed in linger!
  }

  void set_linger(float l_time, float l_grow) {
    linger = l_time;
    linger_grow = l_grow;
  }

  // set target vy, calculate vy_a to achieve it
  // before linger?
  void set_target_vy(float ty) {
    vy_a = (ty - vy) / max_age;
  }
  
  void set_sizes(float min_s, float max_s, float gr, float sr, float sg, float sb, float dd, float tf, float tt) {
    min_size = min_s;
    max_size = max_s;
    r = min_size;
    growth = gr; //growth in m/s
    //
    // if max_size == min_size ... ?
    max_age = (max_size-min_size)/growth;
    /*
    std::string sss = "min_size, max_size, max_age="+std::to_string(min_size)+","+std::to_string(max_size)+","+std::to_string(max_age)+"\n";
    XPLMDebugString( sss.c_str() );
    */
    trans_from = tf; // parameter?
    trans_to = tt; // check if zero, for the linger time, we want very small, but not 0!
    tp = trans_from; // this was missing, have we always started from 0.3?
    //trans_growth = growth/(max_size-min_size)*(trans_from-trans_to);
    trans_growth = (trans_from-trans_to) / max_age; // don't count linger
    //
    smoke_rgb[0] = sr;
    smoke_rgb[1] = sg;
    smoke_rgb[2] = sb;
    //
    draw_delay = dd;
    status = 1; // ignore dd now
    //
    // calculations.
    // to double the speed(s) over the max_age, add Vx / max_age
    // or subtract t0 go to 0
    // maybe count linger time also!!
    // and wind
    // or just a fixed accel. until we reach wind speed!
    vx_a = (wind_x - vx) / max_age; // in max_age, then stay at speed in linger! but if wind 0, motionless...
    // ------>  ? vy_a = 4.0 * (vy / (max_age+linger)); // no wind
    /*
    float vy_desired = -1;
    // if down, we need a height probe...
    // or adapt vy_a so it stops at 0 height?
    vy_a = - (vy - vy_desired) / (max_age); // +linger) removed
    */
    vz_a = (wind_z - vz) / max_age; // in max_age, then stay at speed in linger!
    //
    /*
    vx_a = 4.0 * (vx / (max_age+linger));
    vy_a = 4.0 * (vy / (max_age+linger));
    vz_a = 4.0 * (vz / (max_age+linger));
    */
    /*
    std::string sss = "accel = "+std::to_string(vx_a)+","+std::to_string(vy_a)+","+std::to_string(vz_a)+"\n";
    XPLMDebugString( sss.c_str() );
    */
  }

  /*
    Maybe we should wait one second before drawing the smoke, so we
    don't draw in our. We do update position &c in that time.
    draw_delay, in combination with status... ?
  */
  int update(float delta_t) {
    //lg.xplm("vx, vy, vy="+std::to_string(vx)+","+std::to_string(vy)+","+std::to_string(vz)+"\n");
    if ( status == 0 ) {
      return 0;
    }
    if ( status == 1 ) {
      if ( r >= max_size ) {
	status = 8; // linger
	growth = linger_grow;
	trans_growth = 0; // don't change transparency when lingering, keep at minimum.
	if ( tp < trans_to ) {
	  tp = trans_to; //we tend to overshoot a little
	}
	// x,z 0, keep on rising? wind direction? (available in chute?)
	// stay at wind speed, don't touch rise/sink
	vx_a = 0;
	vy_a = 0;
	vz_a = 0;
	// DEBUG
	/*
	smoke_rgb[0] = 0;
	smoke_rgb[1] = 0;
	smoke_rgb[2] = 0;
	*/
      }
    }
    if ( status == 8 ) { // age of smoke determined by growth until max_size
      linger -= delta_t;
      if ( linger <= 0.0 ) {
	status = 0; // die 
	return status;
      }
    }
    // always, no matter status (which is prolly wrong):
    // size
    if ( status > 0 ) {
      r += growth * delta_t;
      if ( r < 0.0 ) {
	status = 0;
	return status;
      }
      
      // transparency
      tp -= trans_growth*delta_t;
      
      // Position change
      x += vx * delta_t;
      y += vy * delta_t;
      z += vz * delta_t;

      // Height test
      // supply terrain heigth when releasing smoke?
      /*
      float probe_h = height(x, y, z);
      //lg.xplm(std::to_string(probe_h)+","+std::to_string(y)+"\n");
      if ( y <= probe_h ) {
	vy_a = 0.0; // stop all movement? or blow over the ground, new status?
	vy = 0.0;
	vx = 0.0;
	vz = 0.0;
	vx_a = 0.0;
	vz_a = 0.0;
	growth = 0; // dangerous? never stop if not linger? so set linger (status 8)
	//lg.xplm("ground\n");
	status = 8;
      }
      */
      // ADD TEST FOR vy IS DOWNWARDS! OTHERWISE WE DONT TRIGGER HERE
      if ( (status != 8) && (vy < 0.0) && (y < ground_zero) ) { // (y-r) ?
	// just stop going down, all else as normal
	vy_a = 0.0; // stop all movement? or blow over the ground, new status?
	vy = 0.0;
      }
      /* this was the old behaviour, stop all movement
      if ( (status != 8) && (y < ground_zero) ) { // (y-r) ?
	vy_a = 0.0; // stop all movement? or blow over the ground, new status?
	vy = 0.0;
	//
	vx = 0.0; // pn_rand()*0.05; //0.0;
	vz = 0.0; // pn_rand()*0.05; //0.0;
	vx_a = 0.0;
	vz_a = 0.0;
	// take linger growth, then we can set 0, or something else if we want to
	// we could even decrease y with the linger growth, so we don't grow in height
	// ( -> vy = -linger_growth (or half, it's a "radius") )
	growth = linger_grow;
	trans_growth = 0; // normally we are at lowest tp when linger starts
	if ( tp < trans_to ) {
	  tp = trans_to; //we tend to overshoot a little
	}
	//vy = -linger_grow;
	//lg.xplm("ground "+std::to_string(y)+"\n");
	status = 8;	
      }
      */
      
      // Height test

      // Velocity change
      vx += vx_a * delta_t;
      vy += vy_a * delta_t;
      vz += vz_a * delta_t;
      
      // total age
      age += delta_t;
      
    } else {
      // draw delay
      draw_delay -= delta_t;
      if ( draw_delay <= 0.0 ) {
	status = -status; // so we can do this for both 1 and 2
      }
    }
    
    return status;
  }
};

/*
  RisingSmokeBol *b = new RisingSmokeBol();
  b->set_sizes(1.0, 10.0, 0.5,    0.6, 0.6, 0.6,   1.0,    0.8, 0.1);
  b->x = dr_plane_lx;
  b->y = dr_plane_ly;
  //b->vy = 1;
  b->z = dr_plane_lz;
  m_particleList.push_back(b);
*/
 class RisingSmokeBol : public CBol {
 public:

   float life;
   
  // Constructor.
  RisingSmokeBol() {
    status = 1;
    age = 0.0;
    life = 20; // here we take life until 0 ((seconds)
    tp = trans_from;
    r = min_size;
    wind_x = 0.0;
    wind_z = 0.0;
    vx = pn_rand();
    vy = p_rand()*0.2; //1.0;
    vz = pn_rand();
  };

  void set_wind(float wx, float wz) {
    wind_x = wx;
    wind_z = wx;
  }
  
  void set_sizes(float min_s, float max_s, float gr, float sr, float sg, float sb, float dd, float tf, float tt) {
    min_size = min_s;
    max_size = max_s;
    r = min_size;
    growth = gr; // (max_size-min_size) / life; // or as supplied, with max_size

    trans_from = tf; 
    trans_to = tt; 
    tp = trans_from;
    trans_growth = (trans_from-trans_to) / life;
    
    smoke_rgb[0] = sr;
    smoke_rgb[1] = sg;
    smoke_rgb[2] = sb;
    //
    draw_delay = dd;
    status = -1; // -1 is 1 with draw delay. We use distance offset when we emit
    //
    // accelerations, no wind speed
    // set max_speed (terminal velocity) and use wind speed?
    //           *0.1 maybe better than 0.2 ... experiment more.
    vx_a = wind_x*0.1;//(4*wind_x / life);
    vy_a = p_rand()*0.9; // stop rising at life==0, then linger x/z coordinates? or fall down. -p_rand() for AG sparay?
    vz_a = wind_z*0.1; //(4*wind_z / life);
  }

  int update(float delta_t) {
    if ( status == 0 ) {
      return 0;
    }
    if ( status == 1 ) {
      life -= delta_t;

      if ( life <= 0.0 ) {
	status = 0;
	return status;
      }
      
      // size change
      r += growth * delta_t;
      if ( r > max_size ) {
	r = max_size;
	growth = 0;
      }

      /* doesn't look so good.
      if ( life < 10 ) {
	vy_a = -1;
      }
      */

      // Position change
      x += vx * delta_t;
      y += vy * delta_t;
      z += vz * delta_t;
      
      // Velocity change
      vx += vx_a * delta_t;
      vy += vy_a * delta_t;
      vz += vz_a * delta_t;
      
      // Add wind already?
    } else if ( status < 0 ) {
      // draw delay, move, but nothing else... not even age
      // (do we even need to move, the object moves...?)
      draw_delay -= delta_t;
      if ( draw_delay <= 0.0 ) {
	status = -status; // so we can do this for both 1 and 2
      }
      return status;
    }
    // tp goes always down
    tp -= trans_growth*delta_t; // from 1 to 0 in ... seconds, 0 is most transparent
    return status;
  }

};

 
std::vector<CBol*> m_particleList;    // particles for this emitter
std::vector<CBol*>::const_iterator cb_i;

// Update function for drawPhase. No textures in this one!
//
static int DrawBol(int inPhase, int inIsBefore, void* inRefcon) {
  //std::string sss = "inPhase="+std::to_string(inPhase)+"\n";
  //XPLMDebugString( sss.c_str() );
  /*
    int                  inEnableFog,    
    int                  inNumberTexUnits,    
    int                  inEnableLighting,        this one 1?
    int                  inEnableAlphaTesting,    
    int                  inEnableAlphaBlending,    
    int                  inEnableDepthTesting,    
    int                  inEnableDepthWriting);    
  */
  //XPLMSetGraphicsState(1, 0, 1, 1, 1, 1, 1);
  XPLMSetGraphicsState(  0, 1, 0, 1, 1, 1, 0 ); //third one now 0: https://developer.x-plane.com/sdk/XPLMGraphics/
  //XPLMSetGraphicsState(  1, 0, 1, 1, 1, 1, 1); // first one 1?a
  //                           0 for brighter
  //XPLMSetGraphicsState(1, 0, 0, 0, 0, 0, 0); // first one 1?
  // XPLMSetGraphicsState(0, 0, 0, 0, 0, 0, 0); // first one 0?
  glPushMatrix();//mine NEEDED?

  
  //////glEnable(GL_CULL_FACE);
  //////glEnable(GL_DEPTH_TEST);  //is like glDepthMask( GL_FALSE ); // ?

  // WAS ON 2019
  glEnable (GL_BLEND);
  glBlendFunc(  GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA); // more normal
  //
  
  //glBlendFunc(GL_SRC_ALPHA, GL_ONE);
  //glBlendFunc(GL_DST_COLOR, GL_SRC_COLOR);	    
  //glBlendFunc(GL_ONE,       GL_ONE_MINUS_SRC_ALPHA); // very bright/light
  
  //glBlendFunc(GL_SRC_ALPHA, GL_SRC_ALPHA);
  ////glBlendFunc(GL_SRC_ALPHA,  GL_CONSTANT_COLOR );
  
    
  /*
  glDisable(GL_TEXTURE_2D);
  glEnable( GL_BLEND );
  glBlendFunc( GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA );
  */

  float size;
  int segs = 4;// parameter?
  int i, j;
  double x, y, z, z1, z2, R, R1, R2;
  z = cos(PI/segs);
  R = sin(PI/segs);

  /*
  float precalc_x[2*segs+1];
  float precalc_y[2*segs+1];
  for (i = 0; i <= 2*segs; i++) {
    //x = R * cos(-i*2.0*PI/(2*segs));
    precalc_x[i] = R * cos(-i*2.0*PI/(2*segs));
    //y = R * sin(-i*2.0*PI/(2*segs));
    precalc_y[i] = R * sin(-i*2.0*PI/(2*segs));
    }
  */

  for (cb_i = m_particleList.begin(); cb_i != m_particleList.end(); ++cb_i) {
    CBol *b = *cb_i;
    if ( b->status <= 0 ) { // <0 for draw_delay
      continue;
    }
    size = b->r;
    glTranslated(b->x, b->y, b->z);
    glColor4f(b->smoke_rgb[0], b->smoke_rgb[1], b->smoke_rgb[2], b->tp);
    
    float r = size;

    /*
    GLUquadricObj* sphere = gluNewQuadric();
    gluQuadricDrawStyle( sphere, GLU_FILL);
    gluQuadricNormals( sphere, GLU_SMOOTH);
    gluQuadricOrientation( sphere, GLU_OUTSIDE);
    //gluQuadricTexture( sphere, GL_TRUE);
    gluSphere( sphere, size, 5, 5);
    */
    
    //Top cap
    glBegin(GL_TRIANGLE_FAN);
    glNormal3d(0,0,1);
    glVertex3d(0,0,r);
    /* moved outside
    z = cos(PI/segs);
    R = sin(PI/segs);
    */
    for (i = 0; i <= 2*segs; i++) {
      x = R * cos(-i*2.0*PI/(2*segs));  // precalc outside loop? 0..7
      y = R * sin(-i*2.0*PI/(2*segs));
      /*
      x = precalc_x[i];
      y = precalc_y[i];
      */
      glNormal3d(x, y, z);
      glVertex3d(r*x, r*y, r*z);
    }
    glEnd();
    
    //Height segments
    for (j = 1; j < segs-1; j++) {
      z1 = cos(j*PI/segs);
      R1 = sin(j*PI/segs);
      z2 = cos((j+1)*PI/segs);
      R2 = sin((j+1)*PI/segs);
      glBegin(GL_TRIANGLE_STRIP);
      for (i = 0; i <= 2*segs; i++) {
	x = R1*cos(-i*2.0*PI/(2*segs));
	y = R1*sin(-i*2.0*PI/(2*segs));
	glNormal3d(x, y, z1);
	glVertex3d(r*x, r*y, r*z1);
	x = R2*cos(-i*2.0*PI/(2*segs));
	y = R2*sin(-i*2.0*PI/(2*segs));
	glNormal3d(x, y, z2);
	glVertex3d(r*x, r*y, r*z2);
      }
      glEnd();
    } 

    //Bottom cap
    glBegin(GL_TRIANGLE_FAN);
    glNormal3d(0,0,-1);
    glVertex3d(0,0,-r);
    /* moved to above
    z = -cos(PI/segs);
    R = sin(PI/segs);
    */
    for (i = 2*segs; i >= 0; i--) {
      x = R*cos(-i*2.0*PI/(2*segs));
      y = R*sin(-i*2.0*PI/(2*segs));
      glNormal3d(x, y, -z);
      glVertex3d(r*x, r*y, -r*z);
    }
    glEnd();
    
    glTranslated(-b->x, -b->y, -b->z); // mine
  }//num

  glDisable(GL_CULL_FACE);
  glDisable(GL_DEPTH_TEST);
  glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);// reset the original one

  glPopMatrix();
  
  glFlush();

  XPLMSetGraphicsState(0, 0, 0, 0, 0, 0, 0); //needed to get texture
  return TRUE;
}
 
// Update function for drawPhase. Textures in this one!
// (http://forums.x-plane.org/index.php?/forums/topic/90039-disable-custom-rendering-from-casting-shadows/)
//
static int DrawBolBillboard(int inPhase, int inIsBefore, void* inRefcon) {
  /*
    This is the prototype for a low level drawing callback. You are
    passed in the phase and whether it is before or after. If you are
    before the phase, return 1 to let x-plane draw or 0 to suppress
    x-plane drawing. If you are after the phase the return value is
    ignored.
    
    Refcon is a unique value that you specify when registering the
    callback, allowing you to slip a pointer to your own data to the
    callback.
    
    Upon entry the OpenGL context will be correctly set up for you and
    OpenGL will be in 'local' coordinates for 3d drawing and panel
    coordinates for 2d drawing. The OpenGL state (texturing, etc.)
    will be unknown.
  */
  if ( inIsBefore ) {
    return 1;
  }
  
  if ( plane_render_type && world_render_type ) {
    int prt = XPLMGetDatai(plane_render_type);
    int wrt = XPLMGetDatai(world_render_type); // seems always 0
    //lg.xplm("prt, wrt="+std::to_string(prt)+", "+std::to_string(wrt)+"\n");
    if ( prt != 2 ) { // we draw only in prt == 2
      return 1; // this pertains to the airplanes phase
    }
    /*
    if ( wrt != 0 ) {
      return 1;
    }
    */
    // weird error when changing HDR mode, when view on object:
    // sh xp64.sh
    // In CG multi-screen space: 853.000000,406.000000
    // xp64.sh: line 9:  3327 Segmentation fault: 11  ./X-Plane.app/Contents/MacOS/X-Plane
    // seem related to object position, which is gone, crashing the cam-view?
  }
  
  /*
    int                  inEnableFog,    
    int                  inNumberTexUnits,        // without this no texture
    int                  inEnableLighting,        // when 0, bright colours, when 1, almost no colour
    int                  inEnableAlphaTesting,    
    int                  inEnableAlphaBlending,    
    int                  inEnableDepthTesting,    // without this draws though plane
    int                  inEnableDepthWriting);    
  */
  glPushMatrix();

  //XPLMSetGraphicsState(0, 1, 0, 0, 1, 1, 0);
  if ( hdr_on ) {
    int hdr = XPLMGetDatai(hdr_on);
    if ( hdr == 1 ) {
      // or other way around?
      //XPLMSetGraphicsState(0, 1, 0, 1, 1, 1, 0);  // hdr on
      XPLMSetGraphicsState(1, 1, 0, 1, 1, 1, 0);  // lightning 1 means black smoke?!
    } else {
      XPLMSetGraphicsState(1, 1, 0, 1, 1, 1, 0);  // hdr off
    }
  } else {
    XPLMSetGraphicsState(1, 1, 0, 1, 1, 1, 0);
  }
  if ( sun_pitch ) {
    float sp = XPLMGetDataf(sun_pitch);
    if ( sp < -4.0 ) {
      XPLMSetGraphicsState(1, 1, 1, 1, 1, 1, 0);
    }
  }
   XPLMSetGraphicsState(
      0,        // No fog, equivalent to glDisable(GL_FOG);
      1,        // One texture, equivalent to glEnable(GL_TEXTURE_2D);
      0,        // No lighting, equivalent to glDisable(GL_LIGHT0);
      0,        // No alpha testing, e.g glDisable(GL_ALPHA_TEST);
      1,        // Use alpha blending, e.g. glEnable(GL_BLEND);
      1,        // No depth read, e.g. glDisable(GL_DEPTH_TEST);
      0);        // No depth write, e.g. glDepthMask(GL_FALSE)
   
  /*
    https://www.opengl.org/discussion_boards/showthread.php/124970-best-blending-mode-for-particles
    glBlendFunc(GL_SRC_COLOR, GL_ONE_MINUS_SRC_COLOR);
    Then do:
    Ratio=(float)RemainingTimeForParticle/LifeSpanOfParticle;
    glColor3f(Ratio, Ratio, Ratio); 
  */
  
  //glEnable(GL_BLEND); // w/o better!
  //glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA); // w/o better!

  //  from http://learnopengl.com/#!Advanced-OpenGL/Blending (no good though)
  //  glBlendFuncSeparate(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA, GL_ONE, GL_ZERO);
    
  //glBlendFunc(GL_SRC_ALPHA, GL_ONE);
  
  //glBlendFunc(GL_ONE, GL_ONE);

  //glBlendFunc(GL_SRC_ALPHA, GL_ONE);
  //glBlendFunc(GL_DST_COLOR, GL_SRC_COLOR);	    
  //glBlendFunc(GL_ONE,       GL_ONE_MINUS_SRC_ALPHA); // very bright/light
  //glBlendFunc(GL_SRC_ALPHA, GL_SRC_ALPHA);
  //glBlendFunc(GL_SRC_ALPHA,  GL_CONSTANT_COLOR );

  float r[3], u[3];
  //camera_directions(r,u,NULL);
  float m[16];
  glGetFloatv(GL_MODELVIEW_MATRIX, m);
  r[0] = m[0];
  r[1] = m[4];
  r[2] = m[8];
  u[0] = m[1];
  u[1] = m[5];
  u[2] = m[9];

  //lg.xplm("cam_dir r="+std::to_string(r[0])+","+std::to_string(r[1])+","+std::to_string(r[2])+"\n");
  //lg.xplm("cam_dir u="+std::to_string(u[0])+","+std::to_string(u[1])+","+std::to_string(u[2])+"\n");

  XPLMBindTexture2d( circle_texture[0], 0 );
  
  glBegin(GL_QUADS); // or GL_QUADS or GL_POINTS, GL_LINES // move outside loop?
  for (cb_i = m_particleList.begin(); cb_i != m_particleList.end(); ++cb_i) {
    CBol *b = *cb_i;
    if ( b->status <= 0 ) { // <0 for draw_delay
      continue;
    }
    //std::string sss = "bol size="+std::to_string(b->r)+"\n";
    //XPLMDebugString( sss.c_str() );

    glColor4f(b->smoke_rgb[0], b->smoke_rgb[1], b->smoke_rgb[2], b->tp);
    //glColor4f(1.0f, 1.0f, 1.0f, 0.4f); // should be ^^^^^

    // -- XPLMBindTexture2d( circle_texture[0], 0 );

    float rf = b->r;
    float uf = b->r;        
    float x = b->x;// + pn_rand()/10.0;
    float y = b->y;
    float z = b->z;// + pn_rand()/10.0;

    // This seems related: http://www.lighthouse3d.com/opengl/billboarding/index.php?billCheat2
    /*
      d  --  c
      |      |
      a  --  b
      a = center - (right + up) * size;
      b = center + (right - up) * size;
      c = center + (right + up) * size;
      d = center - (right - up) * size;
    */
    
    //http://forums.x-plane.org/index.php?/forums/topic/65670-billboardingimposters-in-opengl-in-x-plane-anyone-ever-pulled-it-off/
    /*
    glBegin(GL_LINE_LOOP); // was dbca
    glVertex3f( x - (r[0]+u[0]), y - (r[1]+u[1]), z - (r[2]+u[2]) ); // a
    glVertex3f( x + (r[0]-u[0]), y + (r[1]-u[1]), z + (r[2]-u[2]) ); // b
    glVertex3f( x + (r[0]+u[0]), y + (r[1]+u[1]), z + (r[2]+u[2]) ); // c
    glVertex3f( x - (r[0]-u[0]), y - (r[1]-u[1]), z - (r[2]-u[2]) ); // d
    glEnd();
    */
    
    // -- glBegin(GL_QUADS); // or GL_QUADS or GL_POINTS, GL_LINES // move outside loop?
    //lg.xplm("xyz="+std::to_string(x)+","+std::to_string(y)+","+std::to_string(z)+"\n");
    // was dbca
    // abcd works also, but upside down?
    // we need -0.5,.5 in texture coordinates?
    
    glTexCoord2f(0.0, 0.0); //glVertex3f(x,   y,   z);
    glVertex3f( x - (r[0] + u[0])*uf, y - (r[1] + u[1])*uf, z - (r[2] + u[2])*uf ); // a
    
    glTexCoord2f(1.0, 0.0); //glVertex3f(x+1, y,   z);
    glVertex3f( x - (r[0] - u[0])*uf, y - (r[1] - u[1])*uf, z - (r[2] - u[2])*uf ); // d

    glTexCoord2f(1.0, 1.0); //glVertex3f(x+1, y+1, z);
    glVertex3f( x + (r[0] + u[0])*uf, y + (r[1] + u[1])*uf, z + (r[2] + u[2])*uf ); // c
    
    glTexCoord2f(0.0, 1.0); //glVertex3f(x,   y+1, z);
    glVertex3f( x + (r[0] - u[0])*uf, y + (r[1] - u[1])*uf, z + (r[2] - u[2])*uf ); // b
    
    // -- glEnd();

    // texturing (python) http://forums.x-plane.org/index.php?/forums/topic/78151-strange-opengl-alpha-issue-xp1030b5/
    
  }//num
  glEnd();
  //glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);// reset the original one

  glPopMatrix();
  
  /*
    glFlush();
    
    XPLMSetGraphicsState(0, 0, 0, 0, 0, 0, 0);
  */
  return TRUE;
}

#define GLM_FORCE_DEPTH_ZERO_TO_ONE
static int DrawBolBillboardVK(int inPhase, int inIsBefore, void* inRefcon) {

  if ( inIsBefore ) {
    return 1; // Let XP draw.
  }
    
  /*
    int                  inEnableFog,    
    int                  inNumberTexUnits,        // without this no texture
    int                  inEnableLighting,        // when 0, bright colours, when 1, almost no colour
    int                  inEnableAlphaTesting,    
    int                  inEnableAlphaBlending,    
    int                  inEnableDepthTesting,    // without this draws though plane
    int                  inEnableDepthWriting);    
  */
  glPushMatrix();

  //XPLMSetGraphicsState(0, 1, 0, 0, 1, 1, 0);
  ////XPLMSetGraphicsState(1, 1, 0, 1, 1, 1, 0);
  XPLMSetGraphicsState(
		       0,        // No fog, equivalent to glDisable(GL_FOG);
		       1,        // One texture, equivalent to glEnable(GL_TEXTURE_2D);
		       0,        // No lighting, equivalent to glDisable(GL_LIGHT0);
		       0,        // No alpha testing, e.g glDisable(GL_ALPHA_TEST);
		       1,        // Use alpha blending, e.g. glEnable(GL_BLEND);
		       0,        // No depth read, e.g. glDisable(GL_DEPTH_TEST);
		       0);        // No depth write, e.g. glDepthMask(GL_FALSE)
  //XPLMSetGraphicsState(0, 1, 0, 1, 1, 0, 0); // imgui
  glDisable(GL_CULL_FACE); // imgui
  glEnable(GL_SCISSOR_TEST); // imgui
  
  /* //imgui
    glPushClientAttrib(GL_CLIENT_ALL_ATTRIB_BITS);
    glPushAttrib(GL_ENABLE_BIT | GL_COLOR_BUFFER_BIT | GL_TRANSFORM_BIT);
    glDisable(GL_CULL_FACE);
    glEnable(GL_SCISSOR_TEST);
    glEnableClientState(GL_VERTEX_ARRAY);
    glEnableClientState(GL_TEXTURE_COORD_ARRAY);
    glEnableClientState(GL_COLOR_ARRAY);
    glEnable(GL_TEXTURE_2D);
    
    glMatrixMode(GL_PROJECTION);
    glPushMatrix();
    glScalef(1.0f, -1.0f, 1.0f);
  */

  float r[3], u[3];
  //camera_directions(r,u,NULL);
  float m[16];

  // Get the current modelview matrix, viewport, and projection matrix from X-Plane
  XPLMGetDatavf(gModelviewMatrixRef, m, 0, 16); //imgui

  r[0] = m[0];
  r[1] = m[4];
  r[2] = m[8];
  u[0] = m[1];
  u[1] = m[5];
  u[2] = m[9];

  //lg.xplm("cam_dir r="+std::to_string(r[0])+","+std::to_string(r[1])+","+std::to_string(r[2])+"\n");
  //lg.xplm("cam_dir u="+std::to_string(u[0])+","+std::to_string(u[1])+","+std::to_string(u[2])+"\n");

  XPLMBindTexture2d( circle_texture[0], 0 );
  
  glBegin(GL_QUADS); // or GL_QUADS or GL_POINTS, GL_LINES // move outside loop?
  for (cb_i = m_particleList.begin(); cb_i != m_particleList.end(); ++cb_i) {
    CBol *b = *cb_i;
    if ( b->status <= 0 ) { // <0 for draw_delay
      continue;
    }
    //std::string sss = "bol size="+std::to_string(b->r)+"\n";
    //XPLMDebugString( sss.c_str() );

    glColor4f(b->smoke_rgb[0], b->smoke_rgb[1], b->smoke_rgb[2], b->tp);
    //glColor4f(1.0f, 1.0f, 1.0f, 0.4f); // should be ^^^^^

    // -- XPLMBindTexture2d( circle_texture[0], 0 );

    float rf = b->r;
    float uf = b->r;        
    float x = b->x;// + pn_rand()/10.0;
    float y = b->y;
    float z = b->z;// + pn_rand()/10.0;
    //z = z * 2 -1; 
    
    // This seems related: http://www.lighthouse3d.com/opengl/billboarding/index.php?billCheat2
    /*
      d  --  c
      |      |
      a  --  b
      a = center - (right + up) * size;
      b = center + (right - up) * size;
      c = center + (right + up) * size;
      d = center - (right - up) * size;
    */
    
    //http://forums.x-plane.org/index.php?/forums/topic/65670-billboardingimposters-in-opengl-in-x-plane-anyone-ever-pulled-it-off/
    /*
    glBegin(GL_LINE_LOOP); // was dbca
    glVertex3f( x - (r[0]+u[0]), y - (r[1]+u[1]), z - (r[2]+u[2]) ); // a
    glVertex3f( x + (r[0]-u[0]), y + (r[1]-u[1]), z + (r[2]-u[2]) ); // b
    glVertex3f( x + (r[0]+u[0]), y + (r[1]+u[1]), z + (r[2]+u[2]) ); // c
    glVertex3f( x - (r[0]-u[0]), y - (r[1]-u[1]), z - (r[2]-u[2]) ); // d
    glEnd();
    */
    
    // -- glBegin(GL_QUADS); // or GL_QUADS or GL_POINTS, GL_LINES // move outside loop?
    //lg.xplm("xyz="+std::to_string(x)+","+std::to_string(y)+","+std::to_string(z)+"\n");
    // was dbca
    // abcd works also, but upside down?
    // we need -0.5,.5 in texture coordinates?
    
    glTexCoord2f(0.0, 0.0); //glVertex3f(x,   y,   z);
    float zc = z - (r[2] + u[2])*uf;
    //zc = 2 * zc - 1;
    glVertex3f( x - (r[0] + u[0])*uf, y - (r[1] + u[1])*uf, zc ); // a

    zc = z - (r[2] - u[2])*uf;
    //zc = 2 * zc - 1;
    glTexCoord2f(1.0, 0.0); //glVertex3f(x+1, y,   z);
    glVertex3f( x - (r[0] - u[0])*uf, y - (r[1] - u[1])*uf, zc ); // d

    zc = z + (r[2] + u[2])*uf;
    //zc = 2 * zc - 1;
    glTexCoord2f(1.0, 1.0); //glVertex3f(x+1, y+1, z);
    glVertex3f( x + (r[0] + u[0])*uf, y + (r[1] + u[1])*uf, zc ); // c

    zc = z + (r[2] - u[2])*uf;
    //zc = 2 * zc - 1;
    glTexCoord2f(0.0, 1.0); //glVertex3f(x,   y+1, z);
    glVertex3f( x + (r[0] - u[0])*uf, y + (r[1] - u[1])*uf, zc ); // b
    
    // -- glEnd();

    // texturing (python) http://forums.x-plane.org/index.php?/forums/topic/78151-strange-opengl-alpha-issue-xp1030b5/
    
  }//num
glEnd();
//glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);// reset the original one

  glPopMatrix();

  /*
  glFlush();

  XPLMSetGraphicsState(0, 0, 0, 0, 0, 0, 0);
  */
  return TRUE;
}


  
}//namespace
#endif

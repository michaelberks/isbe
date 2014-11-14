#ifndef NCM_DMK_CAMERA_H
#define NCM_DMK_CAMERA_H

#include <nailfold/qtools/ncm_cameras/ncm_camera.h>
#include <nailfold/qtools/ncm_cameras/dmk_listener.h>
#include <nailfold/qtools/ncm_cameras/ncm_dmk_grab_frame_filter.h>
#include "tisudshl.h"
#include "grabber.h"

class ncm_dmk_camera : public ncm_camera
{
	Q_OBJECT
// INTERFACE

public:
  // No member variables here, please

  //: Default ctor
  ncm_dmk_camera(/*T param1 = def1, T param2 = def2, ...*/);

  // The 'Big Three'

  ////: Copy ctor
  //ncm_dmk_camera(ncm_dmk_camera const& rhs);
  ////: dtor
  //~ncm_dmk_camera();
  ////: Assignment operator
  //ncm_dmk_camera& operator=(ncm_dmk_camera const& rhs); // assignment


  //virtual void connect();
  //virtual void disconnect();

	~ncm_dmk_camera();

	virtual bool is_connected() /*= 0*/;
	virtual bool connect_device(/*bool use_existing = true*/) /*= 0*/;
	virtual bool connect_device(bool use_existing) /*= 0*/;
	virtual bool disconnect_device() /*= 0*/;
	virtual bool start_live() /*= 0*/;
	virtual bool stop_live() /*= 0*/;
	virtual bool is_live() /*= 0*/;

	virtual vcl_string get_camera_type() /*= 0*/;
	virtual vcl_string get_camera_id() /*= 0*/;

  virtual void capture_image(vil_image_view<vxl_byte>& dest) /*= 0*/;

	virtual bool get_available_FPS(vcl_vector<double>& available_FPS) /*= 0*/;
	virtual bool get_exposure_range(double& min_exposure, double& max_exposure) /*= 0*/;
	virtual bool get_gain_range(double& min_gain, double& max_gain) /*= 0*/;

	virtual bool get_current_FPS(double& curr_FPS) /*= 0*/;
	virtual bool get_current_exposure(double& curr_exposure) /*= 0*/;
	virtual bool get_current_gain(double& curr_gain) /*= 0*/;
	virtual double get_current_FPS() /*= 0*/;
	virtual double get_current_exposure() /*= 0*/;
	virtual double get_current_gain() /*= 0*/;

	virtual bool is_auto_gain() /*= 0*/;
	virtual bool is_auto_exposure() /*= 0*/;

	virtual bool set_FPS(double new_FPS) /*= 0*/;
	virtual bool set_exposure(double new_exposure) /*= 0*/;
	virtual bool set_gain(double new_gain) /*= 0*/;

	virtual bool switch_auto_gain(bool auto_on) /*= 0*/; 
	virtual bool switch_auto_exposure(bool auto_on) /*= 0*/; 

	virtual void set_queue( ncm_video_frame_queue *_queue_ );

  //: Other operators
  //bool operator==(ncm_dmk_camera const& rhs); // equality
  //bool operator!=(ncm_dmk_camera const& rhs); // inequality

  //ncm_dmk_camera* operator->(); // dereference { return p_; }
  //ncm_dmk_camera& operator*(); // dereference { return *p_; }

  //: Functions for binary IO
  //ncm_dmk_camera* clone() const; // virtual copy constructor
  //short version() const;
  //vcl_string is_a() const;
  //bool is_class(const vcl_string& class_name) const;
  //void b_write(vsl_b_ostream& os) const;
  //void b_read(vsl_b_istream& is) const;
  //void print_summary(vcl_ostream& os) const;

signals:

		void frame_ready();


// IMPLEMENTATION

protected:
  // Members and functions visible to objects of this class and derived classes

private:
  // Members and functions visible only to objects of this class
		
	// The instance of the Grabber class.
	DShowLib::Grabber	dmk_grabber_;

	//Handler for images obtained by the grabber
  DShowLib::FrameHandlerSink::tFHSPtr dmk_sink_;

	//Filter to divert frames from DMK sink into our own memory buffer
	ncm_dmk_grab_frame_filter frame_filter_;

	tIVCDPropertyItemsPtr pVCD_props_;


  //: Friend classes
  //friend ...

  //: Static member variables
  //static int member_;
  //static const const_member_ = 42;

  //: Prohibited member functions
  // Put anything here if you want a compile time error when the user tries to
  // call them (e.g. operator+ for two Dates - nonsensical)
  // especially automatically generated functions:
  //   ncm_dmk_camera(); // default constructor
  //   ncm_dmk_camera(ncm_dmk_camera const& rhs); // copy constructor
  //   ~ncm_dmk_camera(); // destructor
  //   ncm_dmk_camera& operator= (const ncm_dmk_camera& rhs); // assignment
  //   ncm_dmk_camera* operator&(); // address of
  //   ncm_dmk_camera const* operator&() const; // address of (const)
};

//: Define inline functions here
//inline
//T inline_fn(params)
//{
//  ...
//}

//inline
//bool ncm_dmk_camera::operator!= (ncm_dmk_camera const& rhs)
//{ return !(*this==rhs); }

#endif // NCM_DMK_CAMERA_H
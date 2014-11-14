#ifndef NCM_LUMENERA_CAMERA_H
#define NCM_LUMENERA_CAMERA_H

#include <nailfold/ncm_camera.h>

class ncm_lumenera_camera : public ncm_camera 
{
// INTERFACE

public:
  // No member variables here, please

  //: Default ctor
  ncm_lumenera_camera(/*T param1 = def1, T param2 = def2, ...*/);

  // The 'Big Three'

  ////: Copy ctor
  //ncm_lumenera_camera(ncm_lumenera_camera const& rhs);
  ////: dtor
  //~ncm_lumenera_camera();
  ////: Assignment operator
  //ncm_lumenera_camera& operator=(ncm_lumenera_camera const& rhs); // assignment


  virtual void connect();
  virtual void disconnect();

  virtual void capture_image(vil_image_view<vxl_byte>& dest);

  //: Other operators
  //bool operator==(ncm_lumenera_camera const& rhs); // equality
  //bool operator!=(ncm_lumenera_camera const& rhs); // inequality

  //ncm_lumenera_camera* operator->(); // dereference { return p_; }
  //ncm_lumenera_camera& operator*(); // dereference { return *p_; }

  //: Functions for binary IO
  //ncm_lumenera_camera* clone() const; // virtual copy constructor
  //short version() const;
  //vcl_string is_a() const;
  //bool is_class(const vcl_string& class_name) const;
  //void b_write(vsl_b_ostream& os) const;
  //void b_read(vsl_b_istream& is) const;
  //void print_summary(vcl_ostream& os) const;



// IMPLEMENTATION

protected:
  // Members and functions visible to objects of this class and derived classes

private:
  // Members and functions visible only to objects of this class
  //: Friend classes
  //friend ...

  //: Static member variables
  //static int member_;
  //static const const_member_ = 42;

  //: Prohibited member functions
  // Put anything here if you want a compile time error when the user tries to
  // call them (e.g. operator+ for two Dates - nonsensical)
  // especially automatically generated functions:
  //   ncm_lumenera_camera(); // default constructor
  //   ncm_lumenera_camera(ncm_lumenera_camera const& rhs); // copy constructor
  //   ~ncm_lumenera_camera(); // destructor
  //   ncm_lumenera_camera& operator= (const ncm_lumenera_camera& rhs); // assignment
  //   ncm_lumenera_camera* operator&(); // address of
  //   ncm_lumenera_camera const* operator&() const; // address of (const)
};

//: Define inline functions here
//inline
//T inline_fn(params)
//{
//  ...
//}

//inline
//bool ncm_lumenera_camera::operator!= (ncm_lumenera_camera const& rhs)
//{ return !(*this==rhs); }

#endif // NCM_LUMENERA_CAMERA_H
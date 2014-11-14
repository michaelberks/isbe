#ifndef __ncm_apex_h__
#define __ncm_apex_h__

//:
// \file
// \brief An apex of a blood vessel (capillary)
// \author Phil Tresadern
//
// \verbatim
//  Modifications
// \endverbatim

#include <vcl_vector.h>
#include <vcl_string.h>

#include <vgl/vgl_point_2d.h>
#include <vgl/vgl_vector_2d.h>

// forward declarations
class vsl_b_ostream;
class vsl_b_istream;
class vsl_ostream;
class vsl_istream;
class ncm_vessel;

class ncm_apex
{
//  INTERFACE

public:
  //  No member variables here, please

  //: Constructors
  ncm_apex(ncm_vessel& parent_vessel,
           double x, double y);

  //  The 'Big Three'

  //: Destructor
  //~ncm_apex();
  
  //: Copy constructor
  //ncm_apex(ncm_apex const& rhs);

  //: Assignment operator
  //ncm_apex& operator=(ncm_apex const& rhs); // assignment

  //: Set coordinates of the inner point of the apex
  void set_inner_point(double x, double y);

  //: Set coordinates of the outer point of the apex
  void set_outer_point(double x, double y);

  //: Coordinates of the apex centre
  vgl_point_2d<double> centre_point() const;

  //: Coordinates of the inner point of the apex
  vgl_point_2d<double> inner_point() const;

  //: Coordinates of the outer point of the apex
  vgl_point_2d<double> outer_point() const;

  //: Width of the vessel at the apex
  double width() const;


  //: Parent vessel
  ncm_vessel& parent_vessel();
  ncm_vessel const& parent_vessel() const;

  //: Uncertainty (be more specific)
  double uncertainty() const;

  //: Check that apex properties are valid (e.g. that it crosses the vessel)
  bool is_valid() const;


  //  Other operators

  //bool operator==(ncm_apex const& rhs); // equality
  //bool operator!=(ncm_apex const& rhs); // inequality

  //  Functions for binary IO
  short version() const;
  vcl_string is_a() const;
  bool is_class(const vcl_string& class_name) const;
  //ncm_apex* clone() const;

  //: Functions for text IO
  void t_write(vcl_ostream& os) const;
  void t_read(vcl_istream& is);
  

  // Friends
  // These functions call protected member functions (e.g. print_summary()) and
  // therefore must be 'friend's.

  //: Binary file stream output operator for class reference
  friend void vsl_b_write(vsl_b_ostream& bfs, const ncm_apex& b);

  //: Binary file stream input operator for class reference
  friend void vsl_b_read(vsl_b_istream& bfs, ncm_apex& b);

  //: Stream output operator for class reference
  friend vcl_ostream& operator<<(vcl_ostream& os,const ncm_apex& b);

  //: Stream output operator for class pointer
  friend vcl_ostream& operator<<(vcl_ostream& os,const ncm_apex* b);



//  IMPLEMENTATION

protected:
  // Members and functions visible to objects of this class and derived classes
  void b_write(vsl_b_ostream& os) const;
  void b_read(vsl_b_istream& is);
  void print_summary(vcl_ostream& os) const;

private:
  // Members and functions visible only to objects of this class
 
  //: Reference to parent vessel object (as the vessel *should* exist)
  ncm_vessel& vessel_;

  //: Offset of apex centre from vessel anchor point
  vgl_vector_2d<double> centre_offset_;

  //: Offset of inner apex point from vessel anchor point
  vgl_vector_2d<double> inner_offset_;

  //: Offset of outer apex point from vessel anchor point
  vgl_vector_2d<double> outer_offset_;
};

#endif // __ncm_apex_h__
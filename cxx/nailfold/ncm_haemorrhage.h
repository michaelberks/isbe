#ifndef __ncm_haemorrhage_h__
#define __ncm_haemorrhage_h__

//:
// \file
// \brief A haemorrhage
// \author Phil Tresadern
//
// \verbatim
//  Modifications
// \endverbatim

#include <vcl_vector.h>
#include <vcl_string.h>

#include <vgl/vgl_point_2d.h>

// forward declarations
class vsl_b_ostream;
class vsl_b_istream;
class vsl_ostream;
class vsl_istream;

class ncm_annotation;
class ncm_vessel;

class ncm_haemorrhage
{
//  INTERFACE

public:
  //  No member variables here, please

  //: Constructors
  ncm_haemorrhage(ncm_annotation& annotation,
                  double x, double y,
                  ncm_vessel* parent_vessel = NULL);

  //  The 'Big Three'

  //: Destructor
  //~ncm_apex();
  
  //: Copy constructor
  //ncm_apex(ncm_apex const& rhs);

  //: Assignment operator
  //ncm_apex& operator=(ncm_apex const& rhs); // assignment

  //: Clear existing data
  void clear();

  //: Parent vessel
  ncm_vessel const* parent_vessel();

  //: Reference point
  void set_anchor(double x, double y);
  const vgl_point_2d<double>& anchor() const;


  //  Outline methods

  //: Add a point to the outline
  void add_outline_point(const vgl_point_2d<double>& point);
  void add_outline_point(double x, double y);

  //: Number of points in outline
  unsigned n_points() const;

  const vgl_point_2d<double>* outline_point(unsigned i) const;

  //: Whether the outline is defined
  bool outline_is_defined() const;

  //: Delete all points from the outline
  void delete_outline();


  //: Attach to a specific vessel (possibly NULL)
  void attach_to_vessel(ncm_vessel const* vessel);


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
  friend void vsl_b_write(vsl_b_ostream& bfs, const ncm_haemorrhage& b);

  //: Binary file stream input operator for class reference
  friend void vsl_b_read(vsl_b_istream& bfs, ncm_haemorrhage& b);

  //: Stream output operator for class reference
  friend vcl_ostream& operator<<(vcl_ostream& os,const ncm_haemorrhage& b);

  //: Stream output operator for class pointer
  friend vcl_ostream& operator<<(vcl_ostream& os,const ncm_haemorrhage* b);



//  IMPLEMENTATION

protected:
  // Members and functions visible to objects of this class and derived classes
  void b_write(vsl_b_ostream& os) const;
  void b_read(vsl_b_istream& is);
  void print_summary(vcl_ostream& os) const;

private:
  // Members and functions visible only to objects of this class
 
  //: Annotation to which this vessel belongs
  ncm_annotation& annotation_;

  //: Pointer to vessel from where the blood leaked (possibly NULL)
  ncm_vessel const* vessel_;

  //: Reference point for haemorrhage
  vgl_point_2d<double> anchor_;

  //: Points defining the closed boundary of the haemorrhage
  vcl_vector< vgl_point_2d<double> > outline_points_;
};

#endif // __ncm_haemorrhage_h__
#ifndef __ncm_vessel_h__
#define __ncm_vessel_h__

//:
// \file
// \brief A blood vessel (capillary). This is defined by an anchor (reference)
//        point, a vessel path, one or more apices, and vessel properties such
//        as 'enlarged' or 'giant'. A list of vessels is contained within the
//        ncm_annotation class.
//
//        Vessels should be handles by an ncm_annotation object to ensure
//        that all vessels are properly accounted for. This is implemented by
//        making all ncm_vessel constructors private such that they can only be
//        called by ncm_vessel itself or by the ncm_vessel_handler class (a
//        friend of ncm_vessel). ncm_vessel_handler itself has a function,
//        new_vessel(), that creates a new ncm_vessel and returns a pointer.
//        This function is also private to ensure it can only be called by 
//        either the ncm_vessel_handler itself or by an ncm_annotation (a
//        friend of ncm_vessel_handler). Making the function static means that 
//        you don't need to create an instance of the ncm_vessel_handler at all.
//
//        The same applies for deleting vessels - only the ncm_annotation class
//        can do this.
//
//        This setup ensures that only the ncm_annotation can create new vessels
//        without giving too much access (as would be the case if ncm_annotation
//        were made a friend of ncm_vessel directly). It's pretty restrictive
//        but is likely to reduce mistakes in the long run.
//
// \author Phil Tresadern
//
// \verbatim
//  Modifications
// \endverbatim

#include <vcl_iosfwd.h>
#include <vcl_vector.h>
#include <vcl_string.h>

#include <vgl/vgl_point_2d.h>

#include "ncm_vessel_properties.h"

// forward declarations
class ncm_annotation;
class ncm_apex;
class ncm_vessel_handler;

//: The actual vessel class
class ncm_vessel
{
//  INTERFACE

  // allow ncm_vessel_handler to create and destroy an ncm_vessel
  friend ncm_vessel_handler;

public:
  //  No member variables here, please

  //: Clear the existing data
  void clear();

  //: Flag annotation as modified
  void set_modified();

  //: Get index of this vessel in the associated annotation
  int index() const;

  //: Set the position of the anchor point
  void set_anchor_to(const vgl_point_2d<double>& point);

  //: Get the anchor position
  const vgl_point_2d<double>& anchor() const;

  //: Vessel properties
  ncm_vessel_properties& properties();
  const ncm_vessel_properties& properties() const;

  //: Add a point to the arterial limb
  void add_arterial_point(const vgl_point_2d<double>& point, double width = 0.0);
  void add_arterial_point(double x, double y, double width = 0.0);

  //: Add a point to the venous limb
  void add_venous_point(const vgl_point_2d<double>& point, double width = 0.0);
  void add_venous_point(double x, double y, double width = 0.0);

  //: Whether path is defined
  bool path_is_defined() const;

  //: Whether the path is 'suspicious' (i.e. the vessel anchor to which it is
  //  closest is one other than itself).
  bool path_is_suspicious() const;

  //: Delete all points from the vessel path
  void delete_path();

  //: Trim off the venous end of the vessel from the point closest to (x,y)
  void trim_at(double x, double y);

  //: Number of points on the path
  unsigned n_points() const;

  //: The points along the centreline
  const vgl_point_2d<double>* point(unsigned i) const;

  //: Return pointer to endpoint of vessel for each limb 
  //  Returns NULL if the path is not defined
  const vgl_point_2d<double>* arterial_endpoint() const;
  const vgl_point_2d<double>* venous_endpoint() const;

  //: Return the index of the point nearest to (x,y)
  //  Returns -1 if the path is not defined
  int point_index_nearest_to(double x, double y) const;

  //: Return the point nearest to (x,y)
  //  Returns NULL if the path is not defined
  const vgl_point_2d<double>* point_nearest_to(double x, double y) const;


  //: Normal to the vessel at a given point
  vgl_vector_2d<double> normal_at(double x, double y) const;


  //: The width of the vessel at every point
  const vcl_vector<double>& widths() const;


  //  Apex methods

  //: Add new apex
  ncm_apex* add_apex_at(double x, double y); // const?

  //: Delete apex
  void delete_apex(unsigned i);
  void delete_apex(ncm_apex* const apex_address);

  //: How many apices
  unsigned n_apices() const;

  //: The i'th apex
  //  Returns NULL is i out of range
  const ncm_apex* apex(unsigned i) const;
  ncm_apex* apex(unsigned i);


  //: Set/get minimum/maximum distance between points
  //  Static so that they can be called without a specific instance of a vessel
  static void set_minimum_inter_point_distance(double distance);
  static double minimum_inter_point_distance();
  static void set_maximum_inter_point_distance(double distance);
  static double maximum_inter_point_distance();


  //  Other operators

  bool operator==(ncm_vessel const& rhs); // equality
  bool operator!=(ncm_vessel const& rhs); // inequality

  //  Functions for binary IO
  short version() const;
  vcl_string is_a() const;
  bool is_class(const vcl_string& class_name) const;
  //ncm_vessel* clone() const;

  //: Functions for text IO
  void t_write(vcl_ostream& os) const;
  void t_read(vcl_istream& is);
  

  // Friends
  // These functions call protected member functions (e.g. print_summary()) and
  // therefore must be 'friend's.

  //: Binary file stream output operator for class reference
  friend void vsl_b_write(vsl_b_ostream& bfs, const ncm_vessel& b);

  //: Binary file stream input operator for class reference
  friend void vsl_b_read(vsl_b_istream& bfs, ncm_vessel& b);

  //: Stream output operator for class reference
  friend vcl_ostream& operator<<(vcl_ostream& os,const ncm_vessel& b);

  //: Stream output operator for class pointer
  friend vcl_ostream& operator<<(vcl_ostream& os, const ncm_vessel* b);


//  IMPLEMENTATION

protected:
  // Members and functions visible to objects of this class and derived classes
  void b_write(vsl_b_ostream& os) const;
  void b_read(vsl_b_istream& is);
  void print_summary(vcl_ostream& os) const;

private:
  // Members and functions visible only to objects of this class
 
  //  Until I can see a need for copying, I'll keep these private

  //: Constructors
  ncm_vessel();
  ncm_vessel(const ncm_annotation& parent_annotation, 
             const vgl_point_2d<double>& anchor);
  ncm_vessel(const ncm_annotation& parent_annotation, 
             double anchor_x, double anchor_y);

  //: Copy constructor
  ncm_vessel(ncm_vessel const& rhs);

  //: Assignment operator
  ncm_vessel& operator=(ncm_vessel const& rhs);

  //: Destructor
  ~ncm_vessel();

  int vessel_size() const;
  int vessel_shape() const;

  //: Delete all apices attached to this vessel
  void delete_all_apices();

  //: Annotation to which this vessel belongs
  const ncm_annotation& annotation_;

  //: Reference point for the vessel
  vgl_point_2d<double> anchor_;

  //: Points along the centreline
  vcl_vector< vgl_point_2d<double> > points_;

  //: Vector of widths, one for every point on the centreline
  vcl_vector<double> widths_;

  //: Vector of apices
  vcl_vector<ncm_apex*> apices_;

  //: Vessel description (OR-ed descriptors)
  ncm_vessel_properties properties_;
  
  //  Static member variables that apply to all vessels

  //: Minimum distance between points (default 0)
  static double minimum_inter_point_distance_;

  //: Maximum distance between points (default 100)
  static double maximum_inter_point_distance_;

  //: How many points either side of the nearest point to use when computing
  //  the normal to the vessel (default 1)
  static unsigned normal_smoothing_;


  //: Friend classes
  //friend ...

  //: Static member variables
  //static int member_;
  //static const const_member_ = 42;

  //: Prohibited member functions
  // Put automatically generated functions here if you want a compile time 
  // error when the user tries to call them
  //   ncm_vessel(); // default constructor
  //   ncm_vessel(ncm_vessel const& rhs); // copy constructor
  //   ~ncm_vessel(); // destructor
  //   ncm_vessel& operator= (const ncm_vessel& rhs); // assignment
  //   ncm_vessel* operator&(); // address of
  //   ncm_vessel const* operator&() const; // address of (const)
};

//: Class defined to give ncm_annotation, and only ncm_annotation, the ability
//  to create and destroy instances of the vessel class
class ncm_vessel_handler
{
  // allow ncm_annotation to create and destroy an ncm_vessel
  friend ncm_annotation;

private:

  // call ncm_vessel constructor and return the new 
  static ncm_vessel* new_vessel(const ncm_annotation& parent_annotation, 
                         const vgl_point_2d<double>& anchor);
  
  static ncm_vessel* new_vessel(const ncm_annotation& parent_annotation, 
                         double anchor_x, double anchor_y);

  static void delete_vessel(ncm_vessel* vessel);
};

#endif // __ncm_vessel_h__